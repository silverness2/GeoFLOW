/*
 * profiler.hpp
 *
 *  Created on: May 24, 2019
 *      Author: bflynt
 */

#ifndef PROFILER_PROFILER_HPP_
#define PROFILER_PROFILER_HPP_


#if defined( USE_PROFILE )

#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>



namespace xstd {


class Profiler final {

public:
	using system_clock      = std::chrono::system_clock;
	using steady_clock      = std::chrono::steady_clock;
	using system_duration   = typename system_clock::duration;
	using steady_duration   = typename steady_clock::duration;
	using system_time_point = typename system_clock::time_point;
	using steady_time_point = typename steady_clock::time_point;
	using String          = std::string;
	using LineNumber      = unsigned;
	using SizeType        = std::size_t;
	using Stack           = std::stack<String>;

	class Injection final {
	public:
		Injection() = delete;
		Injection(const Injection&) = delete;
		Injection(Injection&&) = delete;
		Injection& operator=(const Injection&) = delete;
		Injection& operator=(Injection&&) = delete;

		Injection(const String name, const String file, const LineNumber line)
		: system_(system_clock::now()),
		  steady_(steady_clock::now()){
			Profiler::instance().start(name,file,line);
		}

		~Injection(){
			Profiler::instance().stop(
					system_clock::now() - system_,
					steady_clock::now() - steady_);
		}

	private:
		system_time_point system_;
		steady_time_point steady_;
	};

	// Return One and Only Instance of Singleton
	static Profiler& instance(){
		static Profiler instance_;
		return instance_;
	}

	~Profiler(){
		assert( stack_.empty() );
		this->display();
	}

	void start(
			const String     name,
			const String     file,
			const LineNumber line){
		stack_.push(name);
		if( total_.count(name) == 0 ){
			total_.emplace(
					std::piecewise_construct,
					std::forward_as_tuple(name),
					std::forward_as_tuple(name,file,line));
		}
	}

	void stop(
			const system_duration system_dur,
			const steady_duration steady_dur){
		assert( not stack_.empty() );
		const String name = stack_.top();
		stack_.pop();

		// Add to my totals
		auto acc_iter = total_.find(name);
		assert(acc_iter != total_.end());
		(*acc_iter).second.increment_all(system_dur,steady_dur);

		// Subtract my total from my parents exclusive time
		if( not stack_.empty() ){
			const String parent_name = stack_.top();
			auto parent_iter = total_.find(parent_name);
			assert(parent_iter != total_.end());
			(*parent_iter).second.sub_exclusive(system_dur,steady_dur);
		}
	}

	void display(){
		sort_to_disk_();
	}

private:

	class Accumulator final {
	public:
		Accumulator() = delete;
		Accumulator(const Accumulator&) = default;
		Accumulator(Accumulator&&) = default;
		Accumulator& operator=(const Accumulator&) = default;
		Accumulator& operator=(Accumulator&&) = default;

		Accumulator(const String name, const String file, const LineNumber line)
			: name_(name), file_(file), count_(0), line_(line){}

		void increment_all(
				const system_duration& system_dur,
				const steady_duration& steady_dur){
			system_excl_ += system_dur;
			steady_excl_ += steady_dur;
			system_      += system_dur;
			steady_      += steady_dur;
			count_       += 1;
		}

		void sub_exclusive(
				const system_duration& system_dur,
				const steady_duration& steady_dur){
			system_excl_ -= system_dur;
			steady_excl_ -= steady_dur;
		}

		bool operator<(const Accumulator& rhs){
			return (steady_excl_/count_ < rhs.steady_excl_/rhs.count_);
		}

		String          name_;
		String          file_;
		system_duration system_excl_;
		steady_duration steady_excl_;
		system_duration system_;
		steady_duration steady_;
		SizeType        count_;
		LineNumber      line_;
	};


	void sort_to_disk_() const{
		std::vector<Accumulator> vacc;
		vacc.reserve(total_.size());
		for(const auto& acc : total_){
			vacc.push_back(acc.second);
		}
		std::sort(vacc.begin(),vacc.end());
		std::reverse(vacc.begin(),vacc.end());
		const std::string log_filename = "profile.txt";
		std::ofstream log_file;
		log_file.open(log_filename);
		for(const auto& acc : vacc){
			log_file << acc.name_ << "\n";
			log_file << "File: "   << acc.file_ << "\n";
			log_file << "Line: "   << acc.line_ << "\n";
			log_file << "Calls: "  << acc.count_ << "\n";
			log_file << "System: " << acc.system_.count() << "\n";
			log_file << "Steady: " << acc.steady_.count() << "\n";
			log_file << "System Exclusive: " << acc.system_excl_.count() << "\n";
			log_file << "Steady Exclusive: " << acc.steady_excl_.count() << "\n";
			log_file << "\n";
		}
		log_file.close();
	}


	using UnorderedMap = std::unordered_map<String,Accumulator>;

	UnorderedMap total_;
	Stack        stack_;
};

} /* namespace xstd */




/**
 * \def PROFILE(name)
 * \brief A macro to insert a Profiler
 * \param name Name of function to profile
 */
#ifndef PROFILE
//#define CONCAT_INDIRECT(x, y) x ## y
//#define CONCAT(x, y) CONCAT_INDIRECT(x, y)
#define PROFILE(name) xstd::Profiler::Injection local_scope_profiler(name,__FILE__,__LINE__);
#endif

#else  /* #ifndef USE_PROFILE */
/**
 * Macro to ignore Profiler construction
 */
#ifndef PROFILE
#define PROFILE(name) /*--Empty--*/
#endif

#endif /* #ifndef USE_PROFILE */




#endif /* PROFILER_PROFILER_HPP_ */
