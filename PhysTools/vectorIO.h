#ifndef PHYSTOOLS_VECTORIO_H
#define PHYSTOOLS_VECTORIO_H

namespace phys_tools{

///Wraps a vector of doubles for convenient printing
template<typename T>
struct prettyPrint{
private:
	std::vector<T>& v;
public:
	prettyPrint(std::vector<T>& v):v(v){}

	template<typename U>
	friend std::ostream& operator<<(std::ostream&, const prettyPrint<U>&);
	template<typename U>
	friend std::istream& operator>>(std::istream&, prettyPrint<U>&);
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const prettyPrint<T>& p){
	bool first=true;
	os << '{';
	for(auto x : p.v){
		if(!first)
			os << ", ";
		else
			first=false;
		os << x;
	}
	return(os << '}');
}

template<typename T>
std::istream& operator>>(std::istream& is, prettyPrint<T>& p){
	char c;
	is >> c;
	if(c!='{'){
		is.setstate(is.rdstate()|std::ios::failbit);
		return(is);
	}
	bool done=false;
	std::vector<T> buffer;
	while(!done){
		c=is.peek();
		if(is.eof())
			break;
		switch(c){
			case '}':
				is >> c;
				done=true;
				break;
			case ',':
				is >> c;
				//fall through
			default:
			{
				T d;
				is >> d;
				if(is.fail())
					return(is);
				buffer.push_back(d);
			}
				break;
		}
	}
	p.v.swap(buffer);
	return(is);
}

} //namespace phys_tools

#endif
