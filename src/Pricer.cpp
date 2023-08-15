#include "Pricer.h"
#include <cmath>
#include <vector>
// #include <corecrt_math_defines.h>

MonteCarloPricer::MonteCarloPricer(const PayOff& payoff, const double& maturity, const PathSimulator& path,const size_t& number_of_simulations):_path(path.clone()),_payoff(payoff.clone()),_maturity(maturity), _number_of_simulations(number_of_simulations)
{
}

MonteCarloPricer* MonteCarloPricer::clone() const
{
	return new MonteCarloPricer(*this);
}

MonteCarloPricer::MonteCarloPricer(const MonteCarloPricer& pricer):_path(pricer._path->clone()), _number_of_simulations(pricer._number_of_simulations),_payoff(pricer._payoff->clone()),_maturity(pricer._maturity)
{
}

MonteCarloPricer& MonteCarloPricer::operator=(const MonteCarloPricer& pricer)
{
	if (this != &pricer)
	{
		this->_payoff = pricer._payoff->clone();
		this->_maturity = pricer._maturity;
		this->_number_of_simulations = pricer._number_of_simulations;
		this->_path = pricer._path->clone();
	}
	return *this;
}

MonteCarloPricer::~MonteCarloPricer()
{
	delete _path;
	delete _payoff;
}



double MonteCarloPricer::price() const
{
	double s = 0;
	double nan_count = 0;
	for (size_t i = 0; i < _number_of_simulations; i++)
	{
		// s=s+(*_payoff)(_path->path_maturity(_maturity).first)*exp(-_path->get_model()->risk_free_rate() * _maturity);

		double b = exp(-_path->get_model()->risk_free_rate() * _maturity);
		double a = (*_payoff)({_path->path_maturity().first}); 
		s += a * b;
		// 计算每条path的payoff的时候，用花括号将类型转换成vector <double>即可
		// 如果强制类型转换成const vector <double>，会随机生成不同长度而且全是0的vector情况

		// 调试用
		// if (isnan(a)){
		// 	throw "error";
		// 	nan_count++;
		// } else {
		// 	s += a * b;
		// }

		std::cout << "Simulation Path " << i << " : " << a << " * " << b << " = " << s << std::endl;
	}
	return s / (_number_of_simulations - nan_count);
}
