sn_computing: supernode_computing.cpp
	~/intel/bin/icpc -O3 -mavx -mavx2 supernode_computing.cpp Source/lu_gp_sparse_avx2.cpp Source/lu_gp_sparse.cpp -fpermissive -fPIC -lblas -lm -o sn_computing
clean:
	rm -rf  sn_computing
