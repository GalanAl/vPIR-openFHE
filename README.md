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
    
## Test the code
- Test the code:
    - go in the openfhe-development repository
    - mkdir "rep_for_test"
    - cp CMakeLists.User.txt rep_for_test/CmakeLists.txt
    - add the line "add_executable( test_OPENFHE test_OPENFHE.cpp )" in the CMakeLists.txt file
    - 

