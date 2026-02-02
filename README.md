--------------------------------------------------------------------------------
# Test vPIR under openFHE
--------------------------------------------------------------------------------

## OpenFHE install 

- Installation: 
    - git clone https://github.com/openfheorg/openfhe-development.git 
    - mkdir build
    - cd build
    - cmake ..
    - make -j
    - make install (or cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr . && make all install)
    
## Clone this rep

- Go in the openfhe-development repository
- Clone this rep
- cd vPIR-openFHE
- mkdir build
- cd build
- mkdir SizeData
- cmake ..
- make

## Run the program

-./test_OPENFHE C N_q N_a lbd
- C : number of cuts per query ; N_q(=1) : number of ciphertexts queried ; N_a(=1) : number of ciphertexts answered ; lbd(=42) : number of real+verif queries

## Outputs 

- Size of the database wrt the input parameters
- Write the keys, the queries and the answers into SizeData's files
- Estimates the time required for the protocol wrt the cost of each homomorphic basic operation in parallel on the machine
- Give the real time of each part of the protocol (packing, unpacking, evaluation)
- Verifies that the decrypted results corresponds to the expected result

