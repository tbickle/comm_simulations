#include <vector>

//#ifndef BASE_H
//#define BASE_H

namespace gr {
  namespace trellis {

    bool dec2base(unsigned int num, int base, std::vector<int> &s);
    bool dec2bases(unsigned int num, const std::vector<int> &bases, std::vector<int> &s);
    unsigned int base2dec(const std::vector<int> &s, int base);
    unsigned int bases2dec(const std::vector<int> &s, const std::vector<int> &bases);

  } /* namespace trellis */
} /* namespace gr */

//#endif

