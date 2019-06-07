/*
 * unformated_stream.hpp
 *
 *  Created on: Jun 5, 2019
 *      Author: bryan.flynt
 */

#ifndef INCLUDE_XSTD_STREAM_UNFORMATTED_STREAM_HPP_
#define INCLUDE_XSTD_STREAM_UNFORMATTED_STREAM_HPP_

#include <cassert>
#include <fstream>        // std::fstream
#include <string>         // std::string
#include <utility>        // std::get<>(pair)

namespace xstd {

namespace detail {
namespace binary {

	template<typename I, typename T>
	std::size_t read(std::fstream& fs, std::pair<I,T*> pr){
		const auto num_bytes = pr.first * sizeof(T);
		auto data_ptr = reinterpret_cast<char*>(pr.second);
		fs.read(data_ptr, num_bytes);
		return num_bytes;
	}

	template<typename I, typename T, typename ...Args>
	std::size_t read(std::fstream& fs, std::pair<I,T*> pr, Args... args){
		return read(fs,pr) + read(fs,args...);
	}


	template<typename I, typename T>
	std::size_t write(std::fstream& fs, const std::pair<I,T*> pr){
		const auto num_bytes = pr.first * sizeof(T);
		auto data_ptr = reinterpret_cast<const char*>(pr.second);
		fs.write(data_ptr, num_bytes);
		return num_bytes;
	}

	template<typename I, typename T, typename ...Args>
	std::size_t write(std::fstream& fs, const std::pair<I,T*> pr, const Args... args){
		return write(fs,pr) + write(fs,args...);
	}

} // namespace binary

namespace unformatted {

	template<typename I, typename T>
	std::size_t header_size(const std::pair<I,T*> pr){
		return std::get<0>(pr) * sizeof(T);
	}

	template<typename I, typename T, typename ...Args>
	std::size_t header_size(const std::pair<I,T*> pr, const Args... args){
		return header_size(pr) + header_size(args...);
	}

	template<typename Integer>
	Integer read_header(std::fstream& fs){
		Integer value;
		fs.read(reinterpret_cast<char*>(&value), sizeof(Integer));
		return value;
	}

	template<typename Integer>
	void write_header(std::fstream& fs, const Integer& value){
		fs.write(reinterpret_cast<const char*>(&value), sizeof(Integer));
	}

	template<typename ...Args>
	void read(std::fstream& fs, Args&... args){
		using header_type = std::uint32_t;
		auto begin_bytes = read_header<header_type>(fs);
		auto total_bytes = xstd::detail::binary::read(fs,args...);
		auto end_bytes   = read_header<header_type>(fs);
		assert(begin_bytes == total_bytes);
		assert(begin_bytes == end_bytes);
	}

	template<typename ...Args>
	void write(std::fstream& fs, Args&... args){
		using header_type = std::uint32_t;
		auto begin_bytes = header_size(args...);
		write_header<header_type>(fs,begin_bytes);
		auto total_bytes = xstd::detail::binary::write(fs,args...);
		write_header<header_type>(fs,begin_bytes);
		assert(begin_bytes == total_bytes);
	}

} // namespace unformatted
} // namespace detail



class UnformattedStream {

public:

	enum class Status {
		New,
		Old,
		Replace,
		Unknown
	};


	UnformattedStream();
	explicit UnformattedStream(const std::string& filename, Status mode);
	UnformattedStream(const UnformattedStream& other) = delete;
	UnformattedStream(UnformattedStream&& other) = default;

	UnformattedStream& operator=(const UnformattedStream& other) = delete;
	UnformattedStream& operator=(UnformattedStream&& other) = default;


	void open(const std::string& filename, Status mode);

	bool is_open() const;

	void close();

	void swap(UnformattedStream& other);


	template<typename... Args>
	UnformattedStream& read(Args... args);

	template<typename... Args>
	UnformattedStream& write(const Args... args);

private:
	std::fstream  fstream_;
	std::string   filename_;
	Status        mode_;

};

inline
UnformattedStream::UnformattedStream()
	: filename_(""),
	  mode_(Status::Unknown){
}

inline
UnformattedStream::UnformattedStream(const std::string& filename, Status mode)
	: filename_(filename),
	  mode_(mode) {
	this->open(filename,mode);
}

inline
void
UnformattedStream::open(const std::string& filename, Status mode){
	fstream_.exceptions(std::fstream::failbit | std::fstream::badbit);
	filename_ = filename;
	mode_     = mode;
	try {
		switch(mode_) {
		case(Status::New) : {
			//assert(not std::filesystem::exists(filename_));
			fstream_.open(filename_, std::ios::out | std::ios::binary );
			break;
		}
		case(Status::Old) : {
			//assert(std::filesystem::exists(filename_));
			fstream_.open(filename_, std::ios::in | std::ios::out | std::ios::binary );
			break;
		}
		case(Status::Replace) : {
			fstream_.open(filename_, std::ios::trunc | std::ios::out | std::ios::binary );
			break;
		}
		default :
			fstream_.open(filename_, std::ios::in | std::ios::out | std::ios::binary );
		} // end switch
	} // End try
	catch (std::fstream::failure& err) {
        throw "Failed to open file \"%s\": '%s'" , filename_.c_str(), err.what();
    }
}

inline
bool
UnformattedStream::is_open() const{
	return fstream_.is_open();
}

inline
void
UnformattedStream::close(){
	fstream_.close();
}

inline
void
UnformattedStream::swap(UnformattedStream& other){
	using std::swap;
	swap(fstream_,other.fstream_);
	swap(filename_,other.filename_);
	swap(mode_,other.mode_);
}

template<typename... Args>
UnformattedStream&
UnformattedStream::read(Args... args){
	if(fstream_){
		detail::unformatted::read(fstream_, args...);
	}
	else{

	}
	return *this;
}

template<typename... Args>
UnformattedStream&
UnformattedStream::write(const Args... args){
	if(fstream_){
	detail::unformatted::write(fstream_, args...);
	}
	else {

	}
	return *this;
}

} // namespace xstd



#endif /* INCLUDE_XSTD_STREAM_UNFORMATTED_STREAM_HPP_ */
