#ifndef PHYSTOOLS_MCMC_H
#define PHYSTOOLS_MCMC_H

#include <cassert>
#include <cmath>
#include <random>
#include <vector>

#include <PhysTools/optimization/ParameterSet.h>

namespace phys_tools{

struct functionEvaluation{
	std::vector<double> coordinates;
	double value;
};

struct moveProposal{
	std::vector<double> coordinates;
	double logAcceptanceFactor;

	moveProposal()=default;
	moveProposal(std::size_t n):coordinates(n){}
	moveProposal(moveProposal&& mp):
	coordinates(std::move(mp.coordinates)),logAcceptanceFactor(mp.logAcceptanceFactor){}
};

///A Jumper which implements the Metrolpois-Hastings algorithm
///https://doi.org/10.1093%2Fbiomet%2F57.1.97
///\tparam PDist a type which implements sampling and evaluation of a proposal
///              distribution. Must be callable with an `std::vector<double>` to
///              evaluation (a function proportional to) the probability of
///              sampling the point represented by the vector. Must provide a
///              member function `std::vector<double> sample(RNG&, const std::vector<double>&) const`
///              which returns a new set of coodinates sampled from the proposal
///              distribution, optionally taking into account its second argument
///              which represents the 'current point' in the distribution.
template<typename PDist>
struct MetropolisHastings{
	PDist proposal;

	MetropolisHastings(PDist p):proposal(p){}

	template<typename RNG>
	phys_tools::moveProposal operator()(const std::vector<double>& coordinates,
										const std::vector<std::vector<double>>& ensemble,
										RNG& rng) const{
		phys_tools::moveProposal proposedMove;
		proposedMove.coordinates=proposal.sample(rng,coordinates);
		double pProbReverse=proposal(coordinates,proposedMove.coordinates);
		double pProb=proposal(proposedMove.coordinates,coordinates);
		proposedMove.logAcceptanceFactor=log(pProbReverse)-log(pProb);
		return(proposedMove);
	}
};

///A Jumper which implements the 'stretch move' operation from
///Goodman, J. and Weare, J. "Ensemble Samplers with Affine Invariance"
///COMM. APP. MATH. AND COMP. SCI. Vol. 5, No. 1, 2010
struct StretchMove{
	///For a constant a>1, produce random variates with a distribution of
	/// 1/sqrt(x) in the domain [1/a,a]
	struct square_root_distribution{
		double a;
		double norm; //the normalization constant, up to a factor of 2
		mutable std::uniform_real_distribution<double> uni;

		explicit square_root_distribution(double a):
		a(a),norm(1/((sqrt(a)-sqrt(1/a)))),uni(norm*sqrt(1/a),norm*sqrt(a)){
			if(a<=1)
				throw std::domain_error("square_root_distribution requires a>1");
		}

		template<typename RNG>
		double operator()(RNG& rng) const{
			double v=uni(rng)/norm;
			return(v*v);
		}
	};

	square_root_distribution jumpDist;

	explicit StretchMove(double maxScale=2):jumpDist(maxScale){}

	template<typename RNG>
	phys_tools::moveProposal operator()(const std::vector<double>& coordinates,
										const std::vector<std::vector<double>>& ensemble,
										RNG& rng) const{
		assert(!coordinates.empty());
		assert(ensemble.size()>1 && "StretchMove requires a non-trivial ensemble");
		phys_tools::moveProposal proposedMove(coordinates.size());

		//pick a member of the ensemble to move toward, but ensure that it is
		//not the same as the one for which we are making this jump
		std::uniform_int_distribution<std::size_t> idxDist(0,ensemble.size()-1);
		std::size_t idx;
		std::size_t maxTrials=100;
		//small concern: Should the ensmble become entirely degenerate this will
		//become an infinite loop. Deal with this by limiting the number of trials.
		//100 seems safe; with an ensemble of size 2 this has a probability of
		//2^-100 of failing spuriously, and that probability decreases with ensemble
		//size.
		do{
			idx=idxDist(rng);
			if(!--maxTrials)
				throw std::runtime_error("StretchMove failed too many times to find a distinct ensemble member. "
				                         "Ensmeble may have become degenerate.");
		}while(std::equal(coordinates.begin(),coordinates.end(),ensemble[idx].begin()));

		//pick how far to go
		double z=jumpDist(rng);

		//compute where we're proposing to go
		for(std::size_t i=0; i<coordinates.size(); i++)
			proposedMove.coordinates[i]=ensemble[idx][i]
			                            +z*(coordinates[i]-ensemble[idx][i]);
		//and the penalty associated with going there
		proposedMove.logAcceptanceFactor=(coordinates.size()-1)*log(z);

		return(proposedMove);
	}
};

///\param distribution the distribution from which to sample
///\param jumper object which proposes new points to which to move
///\param parameters description of parameter bounds and fixed parameters
///\param nSamples the number of samples to collect
///\param burnIn the number of initial samples to discard
///\param skip the number of samples to discard per sample retained to control
///            autocorrelation
///\param rng the source of random numbers
///\param initialEnsemble initial positions in the sampling space. The number of
///                       entries in this parameter determines the ensemble size.
///\tparam Dist a type which can be called with an `std::vector<double>` to
///             evaluate a distribution at the point represented by the vector
///\tparam Jumper a type which can be called with an `std::vector<double>`
///               representing a current set of coordinates, an
///               `std::vector<std::vector<double>>` representing an ensemble of
///               coordinate sets, and an RNG as a source of randomness
template<typename Dist, typename Jumper, typename RNG>
std::vector<functionEvaluation>
markovSample(const Dist& distribution, const Jumper& jumper,
			 const ParameterSet& parameters,
			 size_t nSamples, size_t burnIn, size_t skip, RNG& rng,
			 const std::vector<std::vector<double>>& initialEnsemble){
	size_t ensembleSize=initialEnsemble.size();
	std::vector<functionEvaluation> results;
	results.reserve(nSamples);
	size_t skipCounter=skip;
	//internally, use coordinates for only the free parameters
	std::vector<std::vector<double>> coordinates(ensembleSize);
	for(size_t i=0; i<ensembleSize; i++){
		coordinates[i].resize(parameters.numberOfFreeParameters());
		parameters.extractFreeParameters(initialEnsemble[i],coordinates[i]);
	}
	std::vector<std::vector<double>> oldCoordinates(ensembleSize);
	std::vector<double> llh(ensembleSize), lastLLH(ensembleSize);
	std::uniform_real_distribution<double> acceptDist(0,1);
	std::vector<double> baseCoordinates=parameters.getParameterValues();
	std::vector<double> proposedCoordinates(parameters.numberOfFreeParameters());

	//collect initial likelihoods
	for(unsigned int i=0; i<ensembleSize; i++){
		std::copy(baseCoordinates.begin(),baseCoordinates.end(),proposedCoordinates.begin());
		parameters.insertFreeParameters(coordinates[i],proposedCoordinates);
		llh[i]=distribution(proposedCoordinates);
	}

	//sample
	while(results.size()<nSamples){
		//store the current ensemble state
		oldCoordinates=coordinates;
		lastLLH=llh;

		//update each member of the ensemble
		for(unsigned int i=0; i<ensembleSize; i++){
			moveProposal proposedMove=jumper(coordinates[i],oldCoordinates,rng);
			//augment internal coordinates with all fixed parameters before
			//evaluating the target distribution
			std::copy(baseCoordinates.begin(),baseCoordinates.end(),proposedCoordinates.begin());
			parameters.insertFreeParameters(proposedMove.coordinates,proposedCoordinates);
			if(!parameters.inBounds(proposedCoordinates))
				continue; //never accept if out of bounds
			double newLLH=distribution(proposedCoordinates);
			double log_ratio=(lastLLH[i]-newLLH)+proposedMove.logAcceptanceFactor;
			bool accept=(log_ratio>=0) || (log_ratio>log(acceptDist(rng)));
			if(accept){
				std::copy(proposedMove.coordinates.begin(),proposedMove.coordinates.end(),coordinates[i].begin());
				llh[i]=newLLH;
			}
		}

		//decide whether to use the current positions
		if(burnIn){
			burnIn--;
			continue;
		}
		if(skipCounter){
			skipCounter--;
			continue;
		}
		else //don't skip!
			skipCounter=skip;

		for(unsigned int i=0; i<ensembleSize && results.size()<nSamples; i++){
			//augment internal coordinates with all fixed parameters before
			//returning to the user
			std::copy(baseCoordinates.begin(),baseCoordinates.end(),proposedCoordinates.begin());
			parameters.insertFreeParameters(coordinates[i],proposedCoordinates);
			results.push_back(functionEvaluation{proposedCoordinates,llh[i]});
		}
	}

	return(results);
}

} //namespace phys_tools

#endif
