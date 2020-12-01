#include "cublas_v2.h"

extern "C"
{
cublasHandle_t cublas_handle();
void cublas_init();
void cublas_destroy();
}
