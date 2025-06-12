--------------------------------------------------------------------------------
# Test vPIR under openFHE
--------------------------------------------------------------------------------

**Authors**:  Daniel S. Roche, Seung Geol Choi, Mayank Varia, Alexis Galan

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
- cmake ..
- make

## Run the program

-./test_OPENFHE X 
- X is the number of cut for the query

