/*
 * command_line.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: bryan.flynt
 */

#ifndef SRC_XSTD_COMMAND_LINE_HPP_
#define SRC_XSTD_COMMAND_LINE_HPP_

#include <map>
#include <string>

namespace xstd {

/// CommandLine parses simple input line commands
/**
 * The class parses the command line arguments into string tokens.
 * It assumes all keys start with a "-" symbol.
 *
 * Keys that start with a single "-" will be treated as multiple
 * key arguments with the following value (if provided) applying
 * to each character in the key.
 *
 * Keys that start with a double "--" will be treated as a single
 * key argument with the following value (if provided) applying
 * only to that full key word.
 *
 * Usage:
 * \code{.cpp}
 * int main(int argc, char* argv[]){
 *
 *    // Parse the command line
 *    CommandLine cmd_line(argc,argv);
 *
 *    // Get the "-i" or "--i" argument (with a default)
 *    std::string file_name = cmd_line.get("i", "default.txt");
 *
 *    // Get the "--verbose" argument (with a default)
 *    std::string verbose = cmd_line.get("verbose", "0");
 *
 *    // Determine if a "z" was provided after a single "-"
 *    // "-z" or "-xyz" or "-xzvf" would all qualify
 *    bool z_flag = cmd_line.exists("z");
 *
 *    return 0;
 * }
 * \endcode
 */
class CommandLine final {

private:
	using key_value_type = std::map<std::string,std::string>;

public:

	/// Empty data structure
	CommandLine() = default;

	/// Copy constructor
	CommandLine(const CommandLine& CL) = default;

	/// Move constructor
	CommandLine(CommandLine&& CL) = default;

	/// Construct with command line arguments
	CommandLine(int argc, char* argv[]);

	/// Destroy the class
	~CommandLine() = default;

	/// Assign the class
	CommandLine& operator=(const CommandLine& CL) = default;

	/// Move the class
	CommandLine& operator=(CommandLine&& CL) = default;

	/// Insert a key-value pair into the data
	void insert(const std::string& key, const std::string& value = "");

	/// Return true if the key exists within the command line arguments
	bool exists(const std::string& key) const;

	/// Returns the value following the search key
	/*
	 * Returns the value that followed the key within the argument list.
	 * If the key existed within the arguments but had no value the method
	 * will return an empty result even if a default is provided.  The default
	 * is returned only when they key did not exist on the command line.
	 */
	std::string get(const std::string& key, const std::string& def = "");


private:
	key_value_type key_value_pairs_;

	// Test if key exists
	bool is_key_(const std::string& str) const;

	// Test if string is a single "-" key
	bool is_short_key_(const std::string& str) const;

	// Strip leading "-" from keys
	std::string get_key_(const std::string& str) const;

	// Generate key-value pairs from a tagged "-" or "--" key
	key_value_type generate_pairs_(const std::string& tagged_key,
			                       const std::string& value) const;
};


inline
CommandLine::CommandLine(int argc, char* argv[]){
	for (int i = 1; i < argc; ++i) {
		if( this->is_key_(std::string(argv[i])) ){
			std::string tagged_key = std::string(argv[i]);

			std::string value;
			if( ((i+1) < argc) && !this->is_key_(std::string(argv[i+1])) ){
				value = std::string(argv[i+1]);
				++i;
			}
			auto kvpairs = generate_pairs_(tagged_key,value);
			key_value_pairs_.insert(kvpairs.begin(),kvpairs.end());
		}
	}
}

inline
void
CommandLine::insert(const std::string& key, const std::string& value){
	key_value_pairs_[key] = value;
}

inline
bool
CommandLine::exists(const std::string& key) const{
	return (key_value_pairs_.count(key) > 0);
}

inline
std::string
CommandLine::get(const std::string& key, const std::string& def) {
	std::string ans(def);
	if( this->exists(key) ){
		ans = key_value_pairs_[key];
	}
	return ans;
}

inline
bool
CommandLine::is_key_(const std::string& str) const {
	return (str.find("-") == 0);
}

inline
bool
CommandLine::is_short_key_(const std::string& str) const {
	return (str.rfind("-") == 0);
}

inline
std::string
CommandLine::get_key_(const std::string& str) const{
	return str.substr(str.rfind("-")+1);
}

inline
typename CommandLine::key_value_type
CommandLine::generate_pairs_(const std::string& tagged_key,const std::string& value)const{
	key_value_type ans;
	auto key = get_key_(tagged_key);
	if( is_short_key_(tagged_key) ){
		for(auto it = key.begin(); it != key.end(); ++it){
			ans[std::string(1,*it)] = value;
		}
	}
	else {
		ans[key] = value;
	}
	return ans;
}


} /* namespace xstd */


#endif /* SRC_XSTD_COMMAND_LINE_HPP_ */
