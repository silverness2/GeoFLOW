#include "cublas_init.h"

extern "C"
{
static cublasHandle_t cublas_handle_;

cublasHandle_t cublas_handle() {
    return cublas_handle_;
}

void cublas_init() {
    cublasCreate(&cublas_handle_);
}

void cublas_destroy() {
    cublasDestroy(cublas_handle_);
}

}
