/*
 * string_to.hpp
 *
 *  Created on: Apr 26, 2019
 *      Author: bflynt
 */

#ifndef STRING_TO_HPP_
#define STRING_TO_HPP_

#include <cstdint>
#include <stdexcept>
#include <string>
#include <limits>



namespace xstd {

template<typename T>
T string_to(const std::string& str, int base = 10);


template<>
bool string_to<bool>(const std::string& str, int base){
  if( str.size() == 0 ){
    throw std::invalid_argument("string_to<bool>");
  }
  bool result = true;
  if( str[0] == 'F' || str[0] == 'f' || str[0] == '0' ){
    result = false;
  }
  return result;
}

template<>
int8_t string_to<int8_t>(const std::string& str, int base){
  const int32_t result = std::stoi(str,nullptr,base);
  if( result > std::numeric_limits<int8_t>::max() ||
	  result < std::numeric_limits<int8_t>::lowest() ){
    throw std::out_of_range("string_to<int8_t>("+str+")");
  }
  return static_cast<int8_t>(result);
}

template<>
int16_t string_to<int16_t>(const std::string& str, int base){
  const int32_t result = std::stoi(str,nullptr,base);
  if( result > std::numeric_limits<int16_t>::max() ||
      result < std::numeric_limits<int16_t>::lowest() ){
    throw std::out_of_range("string_to<int16_t>("+str+")");
  }
  return static_cast<int16_t>(result);
}

template<>
int32_t string_to<int32_t>(const std::string& str, int base){
	return std::stoi(str,nullptr,base);
}

template<>
int64_t string_to<int64_t>(const std::string& str, int base){
	return std::stoll(str,nullptr,base);
}

template<>
uint8_t string_to<uint8_t>(const std::string& str, int base){
	const uint32_t result = std::stoul(str,nullptr,base);
	if( result > std::numeric_limits<uint8_t>::max() ){
		throw std::out_of_range("string_to<uint8_t>("+str+")");
	}
	return static_cast<int8_t>(result);
}

template<>
uint16_t string_to<uint16_t>(const std::string& str, int base){
  const uint32_t result = std::stoul(str,nullptr,base);
	if( result > std::numeric_limits<uint16_t>::max() ){
		throw std::out_of_range("string_to<uint16_t>("+str+")");
	}
	return static_cast<int16_t>(result);
}

template<>
uint32_t string_to<uint32_t>(const std::string& str, int base){
	return std::stoul(str,nullptr,base);
}

template<>
uint64_t string_to<uint64_t>(const std::string& str, int base){
	return std::stoull(str,nullptr,base);
}

template<>
float string_to<float>(const std::string& str, int base){
  return std::stof(str,nullptr);
}

template<>
double string_to<double>(const std::string& str, int base){
  return std::stod(str,nullptr);
}

template<>
long double string_to<long double>(const std::string& str, int base){
  return std::stold(str,nullptr);
}

} /* namespace xstd */

#endif /* STRING_TO_HPP_ */
