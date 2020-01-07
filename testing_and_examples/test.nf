

ref = file(params.reference)

Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { 
    	error "Cannot find any reads matching: ${params.reads}" 
    }  
    .set { read_pairs } 


println ref