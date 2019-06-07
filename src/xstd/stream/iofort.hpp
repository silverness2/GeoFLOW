/*
 * iofort.hpp
 *
 *  Created on: Jun 5, 2019
 *      Author: bflynt
 */

#ifndef INCLUDE_XSTD_STREAM_IOFORT_HPP_
#define INCLUDE_XSTD_STREAM_IOFORT_HPP_

#include <cassert>
#include <fstream>        // std::fstream


namespace xstd {



class iofort {

public:

	using RecordSize = std::int32_t;

	enum class ios {
		read,
		write,
		readwrite
	};

	iofort();
	iofort(const iofort& other) = default;
	iofort(iofort&& other) = default;
	~iofort() = default;
	iofort& operator=(const iofort& other) = default;
	iofort& operator=(iofort&& other) = default;

	iofort(const std::string& filename,
		   const ios mode = ios::readwrite);

	void swap(iofort& other);

	void open();
	void open(const std::string& filename,
			  const ios mode = ios::readwrite);
	bool is_open() const;
	void close();
	void flush();



	template<typename T>
	void read_binary(const std::size_t num_values, T* data_ptr);

	template<typename I, typename T, typename... Pairs>
	void read_binary(std::pair<I,T*> pair);

	template<typename I, typename T, typename... Pairs>
	void read_binary(std::pair<I,T*> pair, Pairs... pairs);

	template<typename T>
	void write_binary(const std::size_t num_values, const T* data_ptr);

	template<typename I, typename T, typename... Pairs>
	void write_binary(const std::pair<I,T*> pair);

	template<typename I, typename T, typename... Pairs>
	void write_binary(const std::pair<I,T*> pair, const Pairs... pairs);


	void skip_record(const std::size_t nrec);

	void goto_record(const std::size_t nrec);

	template<typename T>
	void read_record(const std::size_t num_values, T* data_ptr);

	template<typename I, typename T, typename... Pairs>
	void read_record(std::pair<I,T*> pair);

	template<typename I, typename T, typename... Pairs>
	void read_record(std::pair<I,T*> pair, Pairs... pairs);

	template<typename T>
	void write_record(const std::size_t num_values, const T* data_ptr);

	template<typename I, typename T, typename... Pairs>
	void write_record(const std::pair<I,T*> pair);

	template<typename I, typename T, typename... Pairs>
	void write_record(const std::pair<I,T*> pair, const Pairs... pairs);

	RecordSize get_record_size();

	template<typename I, typename T, typename... Pairs>
	RecordSize calc_record_size(const std::pair<I,T*> pair);

	template<typename I, typename T, typename... Pairs>
	RecordSize calc_record_size(const std::pair<I,T*> pair, const Pairs... pairs);

private:
	std::fstream fstream_;
	std::string  filename_;
	ios          mode_;
};

inline
iofort::iofort()
	: mode_(ios::readwrite){
}

inline
iofort::iofort(const std::string& filename, const ios mode)
	: filename_(filename),
	  mode_(mode){
}

inline
void
iofort::open(){
	assert(filename_.size() > 0);
	switch(mode_) {
	case(ios::read) : {
		fstream_.open(filename_, std::ios::in | std::ios::binary );
		break;
	}
	case(ios::write) : {
		fstream_.open(filename_, std::ios::out | std::ios::binary );
		break;
	}
	default :
		fstream_.open(filename_, std::ios::in | std::ios::out | std::ios::binary );
	} // end switch
	assert( is_open() );
}


inline
void
iofort::open(const std::string& filename, const ios mode){
	filename_ = filename;
	mode_     = mode;
	this->open();
}

inline
bool
iofort::is_open() const{
	return fstream_.is_open();
}

inline
void
iofort::close(){
	fstream_.close();
}

inline
void
iofort::flush(){
	fstream_.flush();
}

template<typename T>
void
iofort::read_binary(const std::size_t num_values, T* data_ptr){
	if( fstream_ ){
		fstream_.read(reinterpret_cast<char*>(data_ptr), num_values*sizeof(T));
	}
	else {
		assert(false);
	}
}

template<typename I, typename T, typename... Pairs>
void
iofort::read_binary(std::pair<I,T*> pair){
	this->read_binary(pair.first,pair.second);
}

template<typename I, typename T, typename... Pairs>
void
iofort::read_binary(std::pair<I,T*> pair, Pairs... pairs){
	this->read_binary(pair.first,pair.second);
	this->read_binary(pairs...);
}



template<typename T>
void
iofort::write_binary(const std::size_t num_values, const T* data_ptr){
	if( fstream_ ){
		fstream_.write(reinterpret_cast<const char*>(data_ptr), num_values*sizeof(T));
	}
	else {
		assert(false);
	}
}

template<typename I, typename T, typename... Pairs>
void
iofort::write_binary(const std::pair<I,T*> pair){
	this->write_binary(pair.first,pair.second);
}

template<typename I, typename T, typename... Pairs>
void
iofort::write_binary(const std::pair<I,T*> pair, const Pairs... pairs){
	this->write_binary(pair.first,pair.second);
	this->write_binary(pairs...);
}

inline
typename iofort::RecordSize
iofort::get_record_size(){
	RecordSize value;
	this->read_binary(1, &value);
	return value;
}

template<typename I, typename T, typename... Pairs>
typename iofort::RecordSize
iofort::calc_record_size(const std::pair<I,T*> pair){
	return static_cast<RecordSize>(pair.first);
}

template<typename I, typename T, typename... Pairs>
typename iofort::RecordSize
iofort::calc_record_size(const std::pair<I,T*> pair, const Pairs... pairs){
	return calc_record_size(pair) + calc_record_size(pairs...);
}

template<typename T>
void
iofort::read_record(const std::size_t num_values, T* data_ptr){
	auto header = this->get_record_size();
	this->read_binary(num_values,data_ptr);
	auto footer = this->get_record_size();
	assert(header==footer);
}

template<typename I, typename T, typename... Pairs>
void
iofort::read_record(std::pair<I,T*> pair){
	this->read_record(pair.first,pair.second);
}

template<typename I, typename T, typename... Pairs>
void
iofort::read_record(std::pair<I,T*> pair, Pairs... pairs){
	auto header = this->get_record_size();
	this->read_binary(pair,pairs...);
	auto footer = this->get_record_size();
	assert(header==footer);
}

template<typename T>
void
iofort::write_record(const std::size_t num_values, const T* data_ptr){
	const RecordSize header = num_values * sizeof(T);
	this->write_binary(1,&header);
	this->write_binary(num_values,data_ptr);
	this->write_binary(1,&header);
}

template<typename I, typename T, typename... Pairs>
void
iofort::write_record(const std::pair<I,T*> pair){
	this->write_record(pair.first,pair.second);
}

template<typename I, typename T, typename... Pairs>
void
iofort::write_record(const std::pair<I,T*> pair, const Pairs... pairs){
	const RecordSize size = calc_record_size(pair,pairs...);
	this->write_binary(1,&size);
	this->write_binary(pair,pairs...);
	this->write_binary(1,&size);
}

inline
void
iofort::skip_record(const std::size_t num_to_skip){
	for(std::size_t i = num_to_skip; i --> 0 ;){
		auto header = this->get_record_size();
		fstream_.seekg(header, std::ios_base::cur);
		auto footer = this->get_record_size();
		assert(header==footer);
	}
}

inline
void
iofort::goto_record(const std::size_t rec_number){
	fstream_.seekg(0, std::ios_base::beg);
	this->skip_record(rec_number);
}


inline
void
iofort::swap(iofort& other){
	using std::swap;
	swap(fstream_,other.fstream_);
	swap(filename_,other.filename_);
	swap(mode_,other.mode_);
}

} /* namespace xstd */





#endif /* INCLUDE_XSTD_STREAM_IOFORT_HPP_ */
