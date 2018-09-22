#ifndef PHYS_TOOLS_PHYSICS_WEIGHTERS_H
#define PHYS_TOOLS_PHYSICS_WEIGHTERS_H

#include <PhysTools/likelihood/weighting.h>

namespace phys_tools{

///Tilt a spectrum by an incremental powerlaw index about a precomputed median energy
///\tparam Event the type of event object to be weighted. Must have a numeric
///              member named 'primaryEnergy'.
///\tparam T the type to use when computing weights
template<typename Event, typename T>
struct powerlawTiltWeighter : public GenericWeighter<powerlawTiltWeighter<Event,T>>{
private:
	double medianEnergy;
	T deltaIndex;
public:
	using result_type=T;
	powerlawTiltWeighter(double me, T dg):
	medianEnergy(me),deltaIndex(dg){}

	result_type operator()(const Event& e) const{
        //only attempt to compute this is we have a primary energy!
        result_type weight=(e.primaryEnergy?pow(e.primaryEnergy/medianEnergy,-deltaIndex):1);
		return(weight);
	}
};

///Apply an exponential cutoff with a given characteristic energy
///\tparam Event the type of event object to be weighted. Must have a numeric
///              member named 'primaryEnergy'.
///\tparam T the type to use when computing weights
template<typename Event, typename T>
struct expCutoffWeighter : public GenericWeighter<expCutoffWeighter<Event,T>>{
private:
	T cutoffEnergy;
public:
	using result_type=T;
	expCutoffWeighter(T ce):
	cutoffEnergy(ce){}

	result_type operator()(const Event& e) const{
        //only attempt to compute this is we have a primary energy!
        result_type weight=(e.primaryEnergy?exp(-e.primaryEnergy/cutoffEnergy):1);
		return(weight);
	}
};

///A log-parabola spectrum
///\tparam Event the type of event object to be weighted. Must have a numeric
///              member named 'primaryEnergy'.
///\tparam T the type to use when computing weights
template<typename Event, typename T>
struct logParabolaWeighter : public phys_tools::GenericWeighter<logParabolaWeighter<Event,T>>{
private:
	double breakEnergy;
	T alpha;
	T beta;
public:
	using result_type=T;
	logParabolaWeighter(double breakEnergy, T alpha, T beta):
		breakEnergy(breakEnergy),alpha(alpha),beta(beta){}

	result_type operator()(const Event& e) const{
		//only attempt to compute this if we have a primary energy!
		result_type weight=pow(e.primaryEnergy/breakEnergy,-alpha -beta*log10(e.primaryEnergy/breakEnergy));
		return(weight);
	}
};

///Compute weights for a spectrum composed of several discrete segments in energy,
///each with its own normalization.
///\tparam Event the type of event object to be weighted. Must have a numeric
///              member named 'primaryEnergy', and an integer member for the index
///              of the segment to which the event belongs
///\tparam T the type to use when computing weights
///\tparam U the type of the index member in the event
template<typename Event, typename T, typename U>
struct segmentedWeighter : public GenericWeighter<segmentedWeighter<Event,T,U>>{
private:
	std::vector<T> segmentWeights;
	U Event::* segmentIndexPtr;
public:
	using result_type=T;
	template<typename Iterator>
	segmentedWeighter(Iterator segmentBegin, Iterator segmentEnd, U Event::* ptr):
	segmentWeights(segmentBegin,segmentEnd),segmentIndexPtr(ptr){}

	result_type operator()(const Event& e) const{
		//only compute if we have weights
		if(segmentWeights.empty())
			return(1.);
		U segmentIndex=e.*segmentIndexPtr;
        //only attempt to compute this is we have a primary energy,
		//give out of range events zero weight
		return(e.primaryEnergy?
			   (segmentIndex>=0 && segmentIndex<segmentWeights.size()?segmentWeights[segmentIndex]:0):
			   1);
	}
};

///Assign a different relative weight to particles and antiparticles
///\tparam Event the type of event object to be weighted. Must have an interger
///              member named 'primaryType' which is positive for particles and
///              negative for anti-particles (following PDG converntion).
///\tparam T the type to use when computing weights
template<typename Event, typename T>
struct antiparticleWeighter : public GenericWeighter<antiparticleWeighter<Event,T>>{
private:
	//1 = equal weighting of particles and antiparticles,
	//0 = zero weight for particles, double weight for antiparticles
	//2 = double weight for particles, zero weight for antiparticles
	T balance;
public:
	using result_type=T;
	antiparticleWeighter(T b):balance(b){}

	result_type operator()(const Event& e) const{
		return((int)e.primaryType<0 ? balance : 2-balance);
	}
};

///Use a member of an event as a weight
///\tparam Event the type of event object to be weighted.
///\tparam T the type to use when computing weights
///\tparam U the type of the Event member to be read as the weight
template<typename Event, typename T, typename U>
struct cachedValueWeighter : public GenericWeighter<cachedValueWeighter<Event,T,U>>{
private:
	U Event::* cachedPtr;
public:
	using result_type=T;
	cachedValueWeighter(U Event::* ptr):cachedPtr(ptr){}
	result_type operator()(const Event& e) const{
		return(result_type(e.*cachedPtr));
	}
};

} //namespace phys_tools

#endif //PHYS_TOOLS_PHYSICS_WEIGHTERS_H
