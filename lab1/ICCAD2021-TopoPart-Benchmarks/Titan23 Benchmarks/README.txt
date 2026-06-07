This folder contains 22 benchmarks.
For each benchmark, there are four parts.
	The first part(line 1) indicates the number of nodes n1;
	The second part(line 2) indicates the number of nets n2;
	The third part(line 3 ~ line 3+n2-1) describes the information of nets:
		Each line in the third part indicates a net.  The corresponding numbers on the line represent the id of nodes in this net.
		Format: <node id> <node id> ...... <node id>
	The fourth part (line 3+n2 ~ last line) describes the information of fixed nodes:
		The i-th line in this part indicates a set of nodes which should be fixed on the FPGA i-1. The corresponding numbers on the line represent the id of nodes in this FPGA.
		Format: <node id> <node id> ...... <node id>
		 
Note: The node index starts from 0.

