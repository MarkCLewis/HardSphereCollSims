#ifndef BIN_INDEX
#define BIN_INDEX

struct BinIndex {
  int i;
	BinIndex operator+(int offset) {
		return BinIndex{i + offset};
	}
	bool valid(int max) { return i >= 0 && i < max; }
};

#endif
