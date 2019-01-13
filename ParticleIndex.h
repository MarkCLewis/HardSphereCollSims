#ifndef PARTICLE_INDEX
#define PARTICLE_INDEX

/**
 * This defines an indexing struct to be used by all populations to
 * refer to particles. This makes indexing particles type-safe instead of
 * of using raw ints. Because C++ has no overhead for wrapping a primitive
 * in a struct, this has no performance implications.
 */

struct ParticleIndex {
	int i;
	bool operator==(ParticleIndex that) { return i == that.i; }
	bool operator!=(ParticleIndex that) { return i != that.i; }
	bool operator>(int a) { return i > a; }
	bool operator>=(int a) { return i >= a; }
	bool operator<(int a) { return i < a; }
	bool operator<=(int a) { return i <= a; }
};

#endif
