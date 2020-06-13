#include "alignment_utils.h"
#include <cstddef>

size_t aligned_alloc_size(size_t count, size_t alignment) {
    if (count <= 0) {
        return 0;
    }
    return ((count - 1)/alignment + 1)*alignment;
}
