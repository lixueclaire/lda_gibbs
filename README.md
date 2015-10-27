LDA Gibbs Sampling
==================================

### parameters
* Count : number of iterations
* alpha : prior parameter
* beta : prior parameter
* K : number of topics

### input
every line in data.txt is a document

### output
output parameters, likelyhood and top 10 words for every topic

### run
		make sampler
		./sampler <data.txt >result
