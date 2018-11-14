#
# Check for C++11x smart pointers
#

check_cxx_source_compiles("
template<class T>
class MyClass{
  T a;
public:
  MyClass(){
  }
  virtual ~MyClass(){
  }
};

// C++11 Support should allow this
template <class T>
using OtherClass = MyClass<T>;

int main(){
  OtherClass<int> C;

  return 0;
}
"
C11_TEMPLATE_TYPEDEF_WORKS ) 
