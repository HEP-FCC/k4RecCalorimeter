#ifndef RECFCCEECALORIMETER_ONNXRUNTIMEUTILITIES_H
#define RECFCCEECALORIMETER_ONNXRUNTIMEUTILITIES_H

#include "onnxruntime_cxx_api.h"

#include <vector>

// convert vector data with given shape into ONNX runtime tensor
template <typename T>
Ort::Value vec_to_tensor(std::vector<T>& data, const std::vector<std::int64_t>& shape,
                         const Ort::MemoryInfo& mem_info) {
  auto tensor = Ort::Value::CreateTensor<T>(mem_info, data.data(), data.size(), shape.data(), shape.size());
  return tensor;
}

#endif
