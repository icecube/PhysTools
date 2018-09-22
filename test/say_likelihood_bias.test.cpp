#include <iostream>
#include <cmath>
#include <random>

#include <PhysTools/histogram.h>
#include <PhysTools/likelihood/likelihood.h>
#include <PhysTools/optimization/lbfgsb/interface.h>

struct Event{
	double x;
	double simNorm;
};

std::ostream& operator<<(std::ostream& os, Event const & e){
	os << e.x;
	os << ' ';
	os << '[' << e.simNorm << ']';
	return(os);
}

struct WeighterMaker{
	const phys_tools::ParameterSet& paramSet;

	template<typename DataType>
	std::function<DataType(const Event&)> operator()(const std::vector<DataType>& params) const{
		assert(params.size()==4);
		DataType backgroundNorm=params[paramSet.getParameterIndex("backgroundNorm")];
		DataType signalNorm=params[paramSet.getParameterIndex("signalNorm")];
		DataType signalPosition=params[paramSet.getParameterIndex("signalPosition")];
		DataType signalWidth=params[paramSet.getParameterIndex("signalWidth")];
		//std::cout << "Visiting point " << backgroundNorm << ',' << signalNorm
		//<< ',' << signalPosition << ',' << signalWidth << std::endl;

		auto weighter=[=](const Event& e)->DataType{
			DataType backgroundWeight=backgroundNorm;
			DataType z=(e.x-signalPosition)/signalWidth;
			DataType signalWeight=20*signalNorm*exp(z*z/-2)/sqrt(8.*atan(1)*signalWidth*signalWidth);
			DataType weight=(backgroundWeight+signalWeight)*e.simNorm;
			return weight;
		};

		return weighter;
	}
};

//no prior
struct TrivialPrior{
	template<typename DataType>
	DataType operator()(const std::vector<DataType>& params) const{
		return(DataType(1));
	}
};

template<typename RNG>
std::vector<Event> generateSimulation(RNG& rng, size_t nEvents, double xMin, double xMax){
	double simNorm=1./nEvents;
	std::uniform_real_distribution<double> dist(xMin,xMax);
	std::vector<Event> events;
	for(size_t i=0; i<nEvents; i++)
		events.push_back(Event{dist(rng),simNorm});
	return(events);
}

template<typename RNG>
std::vector<Event> generateObservation(RNG& rng, size_t nBackground, double xMin, double xMax,
									   size_t nSignal, double signalPosition, double signalWidth){
	std::vector<Event> events;
	nSignal=std::poisson_distribution<size_t>(nSignal)(rng);
	nBackground=std::poisson_distribution<size_t>(nBackground)(rng);
	std::uniform_real_distribution<double> backgroundDist(xMin,xMax);
	for(size_t i=0; i<nBackground; i++)
		events.push_back(Event{backgroundDist(rng),0});
	std::normal_distribution<> signalDist(signalPosition,signalWidth);
	for(size_t i=0; i<nSignal; i++)
		events.push_back(Event{signalDist(rng),0});
	return(events);
}

int main(){
	std::mt19937 rng(52768);

	using BinType=phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>;
	using HistType = phys_tools::histograms::histogram<1,BinType>;
	using phys_tools::histograms::LinearAxis;
	LinearAxis dataAxis(0,.25);

	auto simData=generateSimulation(rng,1e6,0,20);
	auto obsData=generateObservation(rng,1e4,0,20,1e3,12.3,0.6);

	auto bin=[](const std::vector<Event>& data, HistType& hist){
		using phys_tools::histograms::amount;
		for(const auto& e : data)
			hist.add(e.x,amount(std::cref(e)));
	};
	HistType simHist(dataAxis), obsHist(dataAxis);
	bin(simData,simHist);
	bin(obsData,obsHist);

	std::vector<double> seed{1e4,100,1,1};


	phys_tools::ParameterSet params;
	params.addParameter("backgroundNorm");
	params.setParameterValue("backgroundNorm",seed[0]);
	params.setParameterLowerLimit("backgroundNorm",1000);
	params.addParameter("signalNorm");
	params.setParameterValue("signalNorm",seed[1]);
	params.setParameterLowerLimit("signalNorm",0);
	params.addParameter("signalPosition");
	params.setParameterValue("signalPosition",seed[2]);
	params.setParameterLowerLimit("signalPosition", 0);
	params.setParameterUpperLimit("signalPosition",20);
	params.addParameter("signalWidth");
	params.setParameterValue("signalWidth",seed[3]);
	params.setParameterLowerLimit("signalWidth",1e-3);

	WeighterMaker WM{params};

	auto prob=phys_tools::likelihood::makeLikelihoodProblem
	<std::reference_wrapper<const Event>,4>
	 (std::make_tuple(obsHist), {std::make_tuple(simHist)},
	 TrivialPrior(), {1.},
	 phys_tools::likelihood::simpleDataWeighter(),
	 WM,
	 phys_tools::likelihood::SAYLikelihood(),
	 seed,1);
	phys_tools::lbfgsb::LBFGSB_Driver minimizer(params);
	minimizer.setChangeTolerance(1e-6);
	minimizer.setHistorySize(5);
	minimizer.minimize(phys_tools::likelihood::BFGS_Function<decltype(prob)>(prob));
	auto result=minimizer.minimumPosition();
	/*std::cout << "Fit result:";
	for(auto p : result)
		std::cout << ' ' << p;
	std::cout << std::endl;
	std::cout << "LLH: " << minimizer.minimumValue() << std::endl;*/
	if(std::abs(result[0]-1e4)>200)
		std::cout << "Background norm. not well fit: " << result[0] << std::endl;
	if(std::abs(result[1]-1e3)>50)
		std::cout << "Signal norm. not well fit: " << result[1] << std::endl;
	if(std::abs(result[2]-12.3)>.1)
		std::cout << "Signal position not well fit: " << result[2] << std::endl;
	if(std::abs(result[3]-0.6)>.05)
		std::cout << "Signal position not well fit: " << result[3] << std::endl;
}
