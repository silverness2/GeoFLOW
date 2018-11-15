

#include "tbox/pio.hpp"
#include "tbox/io_buffer.hpp"

#include <sstream>

namespace geoflow {
namespace tbox {

//
// Provide initial values to static variables
//
int pio::s_rank = -1;
std::ofstream* pio::s_filestream = NULL;


//
// Define the parallel buffers
//
static IOBuffer pout_buffer;
static IOBuffer perr_buffer;
static IOBuffer plog_buffer;

//
// Define the associated ostream objects
//
std::ostream pio::pout(&pout_buffer);
std::ostream pio::perr(&perr_buffer);
std::ostream pio::plog(&plog_buffer);



void pio::initialize(std::size_t rank){
	s_rank = rank;
	s_filestream = NULL;
  
  // Initialize the standard parallel output stream
  pout_buffer.setActive(s_rank == 0);
  pout_buffer.setPrefixString(std::string());
  pout_buffer.setOutputStream1(&std::cout);
  pout_buffer.setOutputStream2(NULL);
  
  // Initialize the error parallel output stream
  std::stringstream ss;
  ss << "P=" << s_rank << ":";
  std::string buffer = ss.str();
  
  perr_buffer.setActive(true);
  perr_buffer.setPrefixString(buffer);
  perr_buffer.setOutputStream1(&std::cerr);
  perr_buffer.setOutputStream2(NULL);
  
  // Initialize the parallel log file (disabled by default)
  plog_buffer.setActive(false);
  plog_buffer.setPrefixString(std::string());
  plog_buffer.setOutputStream1(NULL);
  plog_buffer.setOutputStream2(NULL);
}


void pio::finalize(){
   std::cout.flush();
   std::cerr.flush();
   shutdownFilestream();
}


void pio::shutdownFilestream(){
   if (s_filestream) {
      s_filestream->flush();
      s_filestream->close();

      delete s_filestream;
      s_filestream = NULL;

      pout_buffer.setOutputStream2(NULL);
      perr_buffer.setOutputStream2(NULL);
      plog_buffer.setOutputStream1(NULL);
      plog_buffer.setActive(false);
   }
}


void pio::logOnlyNodeZero(const std::string& filename){

   // If the filestream was open, then close it and reset streams
   shutdownFilestream();

   // If this is node zero, then open the log stream and redirect output
   if (s_rank == 0) {
      s_filestream = new std::ofstream(filename.c_str());
      if (!(*s_filestream)) {
         delete s_filestream;
         s_filestream = NULL;
         perr << "pio: Could not open log file ``" << filename.c_str() << "''\n";
      } else {
         pout_buffer.setOutputStream2(s_filestream);
         perr_buffer.setOutputStream2(s_filestream);
         plog_buffer.setOutputStream1(s_filestream);
         plog_buffer.setActive(true);
      }
   }
}


void pio::logAllNodes(const std::string& filename){

	// If the filestream was open, then close it and reset streams
   shutdownFilestream();

   // Open the log stream and redirect output
   std::stringstream ss;
   ss << filename << "." << s_rank;
   std::string full_filename = ss.str();

   s_filestream = new std::ofstream(full_filename.c_str());

   if (!(*s_filestream)) {
      delete s_filestream;
      s_filestream = NULL;
      perr << "pio: Could not open log file ``" << full_filename << "''\n";
   } else {
      pout_buffer.setOutputStream2(s_filestream);
      perr_buffer.setOutputStream2(s_filestream);
      plog_buffer.setOutputStream1(s_filestream);
      plog_buffer.setActive(true);
   }
}


void pio::suspendLogging(){
   pout_buffer.setOutputStream2(NULL);
   perr_buffer.setOutputStream2(NULL);
   plog_buffer.setOutputStream1(NULL);
   plog_buffer.setActive(false);
}


void pio::resumeLogging(){
   if (s_filestream) {
      pout_buffer.setOutputStream2(s_filestream);
      perr_buffer.setOutputStream2(s_filestream);
      plog_buffer.setOutputStream1(s_filestream);
      plog_buffer.setActive(true);
   }
}


bool pio::initialized(){
	return (0 <= pio::s_rank);
}

} // namespace tbox
} // namespace geoflow
