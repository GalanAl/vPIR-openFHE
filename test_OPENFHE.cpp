#include <iostream>
#include "openfhe.h"
#include <vector>
#include <algorithm>
#include <random>
#include "omp.h"

using namespace lbcrypto;
using namespace std;

const string DATAFOLDER = "SizeData";

//rotate inplace by k slots to the right (warning: Ctxt is a 2x#slots/2 matrix)
//(0,...,#slot/2-1,#slot/2,...,#slot-1)-->(#slots/2-k,...,#slots/2-k-1,#slots-k,...,#slots-k-1)
void rotate_k(Ciphertext<DCRTPoly>& ctxt, const CryptoContext<DCRTPoly>& cc, const int32_t k)
{
  int32_t sig =k>0?1:-1;
  uint32_t temp = sig*k;
  int32_t count =0;
  while (temp!=0)
  {
    if(temp%2) ctxt = cc->EvalRotate(ctxt, -sig*pow(2,count));
    temp >>= 1;
    count+=1;
  }
  return;
}

//Copy the value in the first slot into the first N slots
void fill_slots(Ciphertext<DCRTPoly>& V,const int32_t N, const CryptoContext<DCRTPoly>& cc)
{
  int32_t r=N;int count =0;long k=1;
  Ciphertext<DCRTPoly> temp2,temp;
  temp = V; bool test = true;
  
  while(r!=0)
  {
    if(r%2==1)
    {
      
      for(;k<pow(2,count);)
      {
        temp2 = temp;
        rotate_k(temp2,cc,k);
        temp+= temp2;  
        k*=2;
      }
      if(test) {V = temp;test = false;}
      else{rotate_k(V,cc,pow(2,count));V+= temp;}
    }
    r>>=1;
    count+=1;
  }
  return;
}

// write index in the base correcponding to the query vector sizes
void query_gen(vector<vector<int64_t>>& query, const int64_t index)
{ 
  int64_t expo = pow(query[0].size(),query.size()); 
  int64_t Q; int64_t R=index;
  for(int i=0;i<int(query.size());++i) 
  {
    expo/=query[i].size();
    Q=R/expo;
    query[i][Q] = int64_t(1); 
    R = R%expo;
  }
  return;
}

// generate random verification vectors
void random_gen(vector<vector<int64_t>>& verif, const CryptoContext<DCRTPoly>& cc)
{
  srand (time(NULL));
  int64_t p = cc->GetCryptoParameters()->GetPlaintextModulus();
  for(int i=0; i<int(verif.size());++i) {for(long j=0;j<int(verif[i].size());++j) verif[i][j] = rand()%p;}  
  return;
}

// pack the query and verification vectors in N_q plaintexts, depending on the choices vector
void packing(vector<Plaintext>& Packs, vector<bool>& choices, vector<vector<int64_t>>& verif, vector<vector<int64_t>>& query, const CryptoContext<DCRTPoly>& cc )
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  int lbd = choices.size(); 
  int S = query[0].size()/Packs.size();
  int32_t row = S*query.size();
  
  for(int i=0;i<int(Packs.size());++i)
  {
    vector<int64_t> int_pack(slots);
    int j =0;
    for(;j<(1+lbd)/2;++j)
    {
      int32_t I=0;
      if(choices[j])
      {
        for(int a=0;a<int(query.size());++a) 
        {
          for(int32_t b=0;b<S;++b) {int_pack[j*row+I]=query[a][i*S+b]; ++I;}
        }
      }
      else
      {
        for(int a=0;a<int(verif.size());++a)
        {
          for(int32_t b=0;b<S;++b) {int_pack[j*row+I]=verif[a][i*S+b]; ++I;}
        }
      }
    }
    int32_t jump = slots/2 -j*row;
    for(;j<lbd;++j)
    {
      int32_t I=0;
      if(choices[j])
      {
        for(int a=0;a<int(query.size());++a) 
        {
          for(int32_t b=0;b<S;++b) {int_pack[j*row+jump+I]=query[a][i*S+b]; ++I;}
        }
      }
      else
      {
        for(int a=0;a<int(verif.size());++a)
        {
          for(int32_t b=0;b<S;++b) {int_pack[j*row+jump+I]=verif[a][i*S+b]; ++I;}
        }
      }
    }
    Packs[i] = cc->MakePackedPlaintext(int_pack);
  }  
  processingTime = TOC(t);
  cout << "Packing time: " << processingTime << "ms" << std::endl;
  return;
}

// extract the coefficients of the num^th vector and copy-paste in row slots
void extract_one(vector<Ciphertext<DCRTPoly>>::iterator beg_V, const vector<Ciphertext<DCRTPoly>>& E_Packs, const int lbd, const int32_t S, const int num, const int32_t row, const CryptoContext<DCRTPoly>& cc)
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  
  #pragma omp parallel for
  for(int32_t i=0;i<S;++i)
  {
    vector<int64_t> mask(slots);
    for (int l=0;l<(1+lbd)/2;++l) {mask[l*row+num*S+i]=1;mask[slots/2+ l*row+i+num*S]=1;}
    for(int k=0;k<int(E_Packs.size());++k)
    {
      *(beg_V+k*S+i) = cc->EvalMult(E_Packs[k] , cc->MakePackedPlaintext(mask));
      rotate_k(*(beg_V+k*S+i),cc,-(i+S*num));
      fill_slots(*(beg_V+k*S+i),row,cc);
    }
  }
  processingTime = TOC(t);
  cout << "Extract_one : " << processingTime/1000 << "s" << std::endl;
  return;
}

// extract the coefficients of each rotations of the last vectors
void extract_one_rot(vector<Ciphertext<DCRTPoly>>::iterator beg_rot, const vector<Ciphertext<DCRTPoly>>& E_Packs, const int lbd, const int32_t S, const int cuts, const CryptoContext<DCRTPoly>& cc)
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  int N_q = E_Packs.size();
  
  #pragma omp parallel for
  for(int32_t i=0;i<S;++i)
  {
    vector<int64_t> mask1(slots),mask2(slots);
    int32_t k=0;
    for(;k<i;++k) 
    {
      for (int l=0;l<(1+lbd)/2;++l)
      {
         mask2[l*S*cuts+(cuts-1)*S+k]=1;mask2[slots/2+ l*S*cuts+(cuts-1)*S+k]=1;
      }
    }
    for(;k<S;++k) 
    {
      for (int l=0;l<(1+lbd)/2;++l)
      {
        mask1[l*S*cuts+(cuts-1)*S+k]=1;mask1[slots/2+ l*S*cuts+(cuts-1)*S+k]=1;
      }
    }
    Ciphertext<DCRTPoly> temp;
    
    for(int k=0;k<N_q;++k)
    {
      temp = cc->EvalMult(E_Packs[(k+cuts-1)%N_q] , cc->MakePackedPlaintext(mask1));
      *(beg_rot+k*S+i) = cc->EvalMult(E_Packs[(k+cuts)%N_q] , cc->MakePackedPlaintext(mask2));
      for(int l=1;l<cuts;++l)
      {
        rotate_k(temp,cc,S);
        temp += cc->EvalMult(E_Packs[(k+cuts-1-l)%N_q] , cc->MakePackedPlaintext(mask1));
        rotate_k(*(beg_rot+k*S+i),cc,S);
        *(beg_rot+k*S+i) += cc->EvalMult(E_Packs[(k+cuts-l)%N_q] , cc->MakePackedPlaintext(mask2));
      }
      rotate_k(temp,cc,-i);
      rotate_k(*(beg_rot+k*S+i),cc,S-i);
      *(beg_rot+k*S+i)+=temp;
      rotate_k(*(beg_rot+k*S+i),cc,-(cuts-1)*S);
    }
  }
  processingTime = TOC(t);
  cout << "Extract_one_rot : " << processingTime/1000 << "s" << std::endl;
  return;
}


void full_extraction(vector<vector<Ciphertext<DCRTPoly>>>& extracts, const vector<Ciphertext<DCRTPoly>>& E_Packs, const int lbd, const int32_t S, const int cuts, const CryptoContext<DCRTPoly>& cc)
{
  for(int i = 0;i<cuts-1;++i)
  {
    vector<Ciphertext<DCRTPoly>> one(E_Packs.size()*S);
    extract_one(one.begin(),E_Packs,lbd, S, i, S*cuts, cc);
    extracts.push_back(one);
  }
  vector<Ciphertext<DCRTPoly>> rot(E_Packs.size()*S);
  extract_one_rot(rot.begin(), E_Packs,lbd, S,cuts,cc);
  extracts.push_back(rot);
  return;
}

// generate a random matrix (arranged by diagonals)
void rand_diag_mat(vector<Plaintext>::iterator beg_M,const int64_t column,const int32_t row, const CryptoContext<DCRTPoly>& cc, const int lbd)
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  
  int64_t p = cc->GetCryptoParameters()->GetPlaintextModulus();
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;

  #pragma omp parallel for
  for(int64_t i = 0;i<column;++i)
  {
      srand (time(NULL));
      vector<int64_t> diag_i(slots);
      
      for(int32_t j = 0;j<row;++j) 
      {
        int64_t r = rand()%p;
        for(int k = 0;k<(lbd+1)/2;++k) {diag_i[k*row+j]=r;diag_i[slots/2+k*row+j]=r;}
      }
      *(beg_M+i) = cc->MakePackedPlaintext(diag_i);
  }
  
  processingTime = TOC(t);
  cout << "matrix_gen: " << processingTime/1000 << "s" << std::endl;
  return;
}

// get the j^th column of a matrix described by its diagonals
void get_col_j(vector<int64_t>::iterator beg_col, vector<Plaintext>::const_iterator beg_M, const int64_t j, const int64_t column, const int64_t row)
{
  Plaintext temp;
  for(long i =0;i<row;++i)
  {
    temp = *(beg_M+((j-i)%column+column)%column);
    *(beg_col+i) = temp->GetPackedValue()[i];
  }
  return;
}
void full_evaluation(vector<Ciphertext<DCRTPoly>>& Answers, vector<vector<Plaintext>>::const_iterator beg_M, vector<vector<Ciphertext<DCRTPoly>>>::const_iterator beg_extract, const int cuts, const int32_t column, const CryptoContext<DCRTPoly>& cc)
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);

  switch(cuts)
  {
    case 1:
    {
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_rot = (*(beg_extract)).begin();
      for(int i=0;i<int(Answers.size());++i)
      {
        vector<Ciphertext<DCRTPoly>> Accu(column);
        vector<Ciphertext<DCRTPoly>>::iterator beg_A = Accu.begin();
        vector<Plaintext>::const_iterator beg_Mi = (*(beg_M+i)).begin();
        #pragma omp parallel for shared(beg_A,beg_Mi)
        for(int32_t A=0;A<column;++A)
        {
                *(beg_A+A) =cc->EvalMult(*(beg_Mi+A),*(beg_rot+A));
        }
        Answers[i] = cc->EvalAddMany(Accu);
      }
      break;
    }
    case 2:
    {
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v = (*(beg_extract)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_rot = (*(beg_extract+1)).begin();
      for(int i=0;i<int(Answers.size());++i)
      {
        vector<Ciphertext<DCRTPoly>> Accu(column);
        vector<Ciphertext<DCRTPoly>>::iterator beg_A = Accu.begin();
        vector<Plaintext>::const_iterator beg_Mi = (*(beg_M+i)).begin();

        #pragma omp parallel for shared(beg_A,beg_Mi,beg_rot,beg_v)
        for(int32_t A=0;A<column;++A)
        {
                *(beg_A+A) = cc->EvalMult(*(beg_Mi),*(beg_rot));
                for(int32_t C=1;C<column;++C)
                {
                        *(beg_A+A) +=cc->EvalMult(*(beg_Mi+C),*(beg_rot+C));
                }
                *(beg_A+A) *= *(beg_v+A);

        }
        Answers[i] = cc->EvalAddMany(Accu);
      }
      break;
    }
    case 3:
    {
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v1 = (*(beg_extract)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v2 = (*(beg_extract+1)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_rot = (*(beg_extract+2)).begin();
      for(int i=0;i<int(Answers.size());++i)
      {
        vector<Ciphertext<DCRTPoly>> Accu(column);
        vector<Ciphertext<DCRTPoly>>::iterator beg_A = Accu.begin();
        vector<Plaintext>::const_iterator beg_Mi = (*(beg_M+i)).begin();
        #pragma omp parallel for shared(beg_A,beg_Mi,beg_rot,beg_v1,beg_v2)
        for(int32_t A=0;A<column;++A)
        {
                vector<Ciphertext<DCRTPoly>> Bccu(column);
                for(int32_t B=0;B<column;++B)
                {
                        Bccu[B] = cc->EvalMult(*(beg_Mi),*(beg_rot));
                        for(int32_t C=1;C<column;++C)
                        {
                                Bccu[B] +=cc->EvalMult(*(beg_Mi+C),*(beg_rot+C));
                        }
                        Bccu[B] *= *(beg_v2+B);
                }
                *(beg_A+A) = cc->EvalMult(*(beg_v1+A),cc->EvalAddMany(Bccu));
        }

        Answers[i] = cc->EvalAddMany(Accu);
      }
      break;
    }
    case 4:
    {
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v1 = (*(beg_extract)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v2 = (*(beg_extract+1)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v3 = (*(beg_extract+2)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_rot = (*(beg_extract+3)).begin();

      for(int i=0;i<int(Answers.size());++i)
      {
        vector<Ciphertext<DCRTPoly>> Accu(column);
        vector<Ciphertext<DCRTPoly>>::iterator beg_A = Accu.begin();
        vector<Plaintext>::const_iterator beg_Mi = (*(beg_M+i)).begin();
        #pragma omp parallel for shared(beg_A,beg_Mi,beg_rot,beg_v1,beg_v2,beg_v3)
        for(int32_t A=0;A<column;++A)
        {
          vector<Ciphertext<DCRTPoly>> Bccu(column);
          for(int32_t B=0;B<column;++B)
          {
            vector<Ciphertext<DCRTPoly>> Cccu(column);
            for(int32_t C=0;C<column;++C)
            {
                Cccu[C] = cc->EvalMult(*(beg_Mi),*(beg_rot));
                for(int32_t D=1;D<column;++D)
                {
                        Cccu[C] +=cc->EvalMult(*(beg_Mi+D),*(beg_rot+D));
                }
                Cccu[C] *= *(beg_v3+C);
            }
            Bccu[B] = cc->EvalMult(cc->EvalAddMany(Cccu),*(beg_v2+B));
          }
          (*(beg_A+A)) = cc->EvalMult(cc->EvalAddMany(Bccu),*(beg_v1+A));
        }
        Answers[i] = cc->EvalAddMany(Accu);
      }
      break;
    }
    case 5:
    {
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v1 = (*(beg_extract)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v2 = (*(beg_extract+1)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v3 = (*(beg_extract+2)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_v4 = (*(beg_extract+3)).begin();
      vector<Ciphertext<DCRTPoly>>::const_iterator beg_rot = (*(beg_extract+4)).begin();

      for(int i=0;i<int(Answers.size());++i)
      {
        vector<Ciphertext<DCRTPoly>> Accu(column);
        vector<Ciphertext<DCRTPoly>>::iterator beg_A = Accu.begin();
        vector<Plaintext>::const_iterator beg_Mi = (*(beg_M+i)).begin();
        #pragma omp parallel for shared(beg_A,beg_Mi,beg_rot,beg_v1,beg_v2,beg_v3,beg_v4)
        for(int32_t A=0;A<column;++A) //TODO AB=0<<column**2 ?
        {
          vector<Ciphertext<DCRTPoly>> Bccu(column);
          for(int32_t B=0;B<column;++B)
          {
            vector<Ciphertext<DCRTPoly>> Cccu(column);
            for(int32_t C=0;C<column;++C)
            {
                vector<Ciphertext<DCRTPoly>> Dccu(column);
                for(int32_t D=0;D<column;++D)
                {
                        Dccu[D] = cc->EvalMult(*(beg_Mi),*(beg_rot));
                        for(int32_t E=1;E<column;++E)
                        {
                                Dccu[D] +=cc->EvalMult(*(beg_Mi+E),*(beg_rot+E));
                        }
                        Dccu[D] *= *(beg_v4+D);
                }
                Cccu[C] = cc->EvalMult(cc->EvalAddMany(Dccu),*(beg_v3+C));;
            }
            Bccu[B] = cc->EvalMult(cc->EvalAddMany(Cccu),*(beg_v2+B));
          }
          (*(beg_A+A)) = cc->EvalMult(cc->EvalAddMany(Bccu),*(beg_v1+A));
        }
        Answers[i] = cc->EvalAddMany(Accu);
      }
      break;
    }
  }

  processingTime = TOC(t);
  cout << "Full_evaluation: " << processingTime/1000 << "s" << std::endl;
  return;
}

bool verification(vector<int64_t>& expected, vector<int64_t>& result, vector<bool>& choice, int32_t sep,const CryptoContext<DCRTPoly>& cc)
{
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  int32_t lbd = choice.size();
  int32_t row = expected.size();
  for(int32_t i =0;i<lbd;++i)
  {
    if(choice[i]==1)
    {
      for(int32_t j=0;j<row;++j)
      {
        if(expected[j]!=result[(i>=lbd/2?slots/2+(i-lbd/2)*sep:i*sep)+j]){cout<<"LOOOSER "<<i<<endl; return false;}
      }
    }
  }
  return true;
}

int main(int argc, char* argv[]) 
{  
    
    if(argc <= 1)
    {
      cerr << "Usage: " << argv[0] << " <number of cuts (1 :: 6)><number of query ctxt (default = 1)><number of answer ctxt (default = 1)><lambda (default = 42)> " << endl;
      exit(0);
    }
    int cuts = (argc>1?atoi(argv[1]):6);
    int N_q = (argc>2?atoi(argv[2]):1);
    int N_a = (argc>3?atoi(argv[3]):1);
    int lbd = (argc>4?atoi(argv[4]):42);
    
    
    
    int64_t numthreads(1);

    #if defined(_OPENMP)
    #pragma omp parallel
    #pragma omp single
    {
      numthreads = omp_get_num_threads();
    }
    #endif

    clog<<"Number of cores : "<<numthreads<<endl;  
      
      
    //********************* SETUP *********************//
    CCParams<CryptoContextBFVRNS> parameters;

    parameters.SetPlaintextModulus(281474978414593);
    parameters.SetMaxRelinSkDeg(3);
    parameters.SetMultiplicativeDepth(5);
    parameters.SetSecurityLevel(HEStd_128_classic);
    

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    
    if (!Serial::SerializeToFile(DATAFOLDER + "/cryptocontext.txt", cc, SerType::BINARY)) {
        std::cerr << "Error writing serialization of the crypto context to "
                     "cryptocontext.txt"
                  << std::endl;
        return 1;
    }
    std::cout << "The cryptocontext has been serialized." << std::endl;
    
    int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    cout<<"nb slot : "<<slots<<endl;
        
    // Initialize the public key containers.
    KeyPair<DCRTPoly> kp = cc->KeyGen();

    std::vector<int32_t> indexList; for(int32_t i=1;i<slots/2;i*=2){ indexList.push_back(int32_t(i));indexList.push_back(int32_t(-i));}
    cc->EvalRotateKeyGen(kp.secretKey, indexList); //Generate all the 2^i rotation keys
    cc->EvalMultKeyGen(kp.secretKey);
    
    
    if (!Serial::SerializeToFile(DATAFOLDER + "/key-public.txt", kp.publicKey, SerType::BINARY)) {
        std::cerr << "Error writing serialization of public key to key-public.txt" << std::endl;
        return 1;
    }
    std::cout << "The public key has been serialized." << std::endl;

    // Serialize the secret key
    if (!Serial::SerializeToFile(DATAFOLDER + "/key-private.txt", kp.secretKey, SerType::BINARY)) {
        std::cerr << "Error writing serialization of private key to key-private.txt" << std::endl;
        return 1;
    }
    std::cout << "The secret key has been serialized." << std::endl;
    
        std::ofstream emkeyfile(DATAFOLDER + "/" + "key-eval-mult.txt", std::ios::out | std::ios::binary);
    if (emkeyfile.is_open()) {
        if (cc->SerializeEvalMultKey(emkeyfile, SerType::BINARY) == false) {
            std::cerr << "Error writing serialization of the eval mult keys to "
                         "key-eval-mult.txt"
                      << std::endl;
            return 1;
        }
        std::cout << "The eval mult keys have been serialized." << std::endl;

        emkeyfile.close();
    }
    else {
        std::cerr << "Error serializing eval mult keys" << std::endl;
        return 1;
    }

    // Serialize the rotation keys
    std::ofstream erkeyfile(DATAFOLDER + "/" + "key-eval-rot.txt", std::ios::out | std::ios::binary);
    if (erkeyfile.is_open()) {
        if (cc->SerializeEvalAutomorphismKey(erkeyfile, SerType::BINARY) == false) {
            std::cerr << "Error writing serialization of the eval rotation keys to "
                         "key-eval-rot.txt"
                      << std::endl;
            return 1;
        }
        std::cout << "The eval rotation keys have been serialized." << std::endl;

        erkeyfile.close();
    }
    else {
        std::cerr << "Error serializing eval rotation keys" << std::endl;
        return 1;
    }
  
    //********************* CLIENT SIDE *********************//
  int32_t S = ((slots/2)/((lbd+1)/2))/cuts; //area per query vector cuts
  int32_t row = S*cuts; //size of ptxt areas = #rows in the matrix
  
  int64_t column = pow(N_q*S,cuts);
  
  cout<<"lbd = "<<lbd<<"; rows = "<<N_a*row<<"; column = "<<column<<"; size small vectors = "<<S<<endl;
  int64_t DB_size = column*N_a*row*6;
  if(DB_size>pow(10,12)){cout<<"max data base size : "<<DB_size/pow(10,12)<<"TB"<<endl;}
  else if (DB_size>pow(10,9)){cout<<"max data base size : "<<DB_size/pow(10,9)<<"GB"<<endl;}
  else if (DB_size>pow(10,6)){cout<<"max data base size : "<<DB_size/pow(10,6)<<"MB"<<endl;}
  else {cout<<"max data base size : "<<DB_size/pow(10,3)<<"KB"<<endl;}
  
  vector<int64_t> length(N_q*S,0);
  vector<vector<int64_t>> query(cuts,length);
  vector<vector<int64_t>> verif(cuts,length);
  
  vector<bool> choice(lbd,true);
  srand (time(NULL));
  
  for(long i=0;i<lbd;++i) choice[i] = bool(rand()%2);
  
  int64_t index = rand()%int64_t(column);
  cout<<"index : "<<index<<" = "<<index%(N_q*S)<<" mod N_q*S"<<endl;
  
  query_gen(query, index);
  random_gen(verif, cc);
  
  Plaintext pack = cc->MakePackedPlaintext(vector<int64_t>(slots,0));
  vector<Plaintext> Packs(N_q,pack);
  Ciphertext<DCRTPoly> E_pack = cc->Encrypt(kp.publicKey, pack);
  vector<Ciphertext<DCRTPoly>> E_Packs(N_q,E_pack);
  
  packing(Packs, choice, verif, query, cc);
    
  for(int i=0;i<N_q;++i) 
  {
    E_Packs[i] = cc->Encrypt(kp.publicKey, Packs[i]); 
    if (!Serial::SerializeToFile(DATAFOLDER + "/query_" +char(i)+ ".txt", E_Packs[i], SerType::BINARY)) 
    {
        std::cerr << "Error writting serialization of query ctxt "<<i<< std::endl;
        return 1;
    }
    std::cout << "The query ciphertexts have been serialized." << std::endl;
  }
  
  vector<int64_t> diag(slots);
  
  // VERIF TIMINGS 
  vector<Ciphertext<DCRTPoly>> accu(numthreads,E_pack);

  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  #pragma omp parallel for
  for (long l=0;l<10000;++l)
  {
    accu[omp_get_thread_num()] += E_Packs[(l+1)%N_q];
  }
  processingTime = TOC(t);
  double hom_add = processingTime/10000;

  Ciphertext<DCRTPoly> Eres = cc->EvalAddMany(accu);
  Plaintext Pres;
  cc->Decrypt(kp.secretKey, Eres, &Pres);
  cout<<bool(Pres==pack);

  accu = vector<Ciphertext<DCRTPoly>>(numthreads,E_pack);

  processingTime=0.0;
  TIC(t);
  #pragma omp parallel for
  for (long l=0;l<10000;++l)
  {
    accu[omp_get_thread_num()] += cc->EvalMult(Packs[l%N_q],E_Packs[l%N_q]);
  }
  processingTime = TOC(t);
  double ptxt_ctxt = processingTime/10000 - hom_add;

  Eres = cc->EvalAddMany(accu);

  cc->Decrypt(kp.secretKey, Eres, &Pres);
  cout<<bool(Pres==pack);

  accu = vector<Ciphertext<DCRTPoly>>(numthreads,E_pack);

  processingTime=0.0;
  TIC(t);
  #pragma omp parallel for
  for (long l=0;l<1000;++l)
  {
    accu[omp_get_thread_num()] += cc->EvalMult(E_Packs[(l+1)%N_q],E_Packs[l%N_q]);
  }
  processingTime = TOC(t);
  double ctxt_ctxt = processingTime/1000-hom_add;
  Eres = cc->EvalAddMany(accu);

  cc->Decrypt(kp.secretKey, Eres, &Pres);
  cout<<bool(Pres==pack)<<endl;

  cout << "ptxt-ctxt product in parallel : " << ptxt_ctxt<< "ms;"<<endl;
  ptxt_ctxt = N_a*double(pow(N_q*S,cuts))*(ptxt_ctxt)/1000;
  cout << "add in parallel : " << hom_add<< "ms;"<<endl;
  double x=0;
  for (int c=0;c<cuts;++c) x+= double(pow(N_q*S,c+1));
  hom_add=N_a*x*(hom_add)/1000;
  cout << "ctxt-ctxt product in parallel : " << ctxt_ctxt<< "ms;"<<endl;
  ctxt_ctxt=N_a*(x/(N_q*S))*ctxt_ctxt/1000;
  cout<<"estimated ptxt-ctxt total time : "<<ptxt_ctxt<<" s" << endl;
  cout<<"estimated add total time : "<<hom_add<<" s" << endl;
  cout<<"estimated ctxt-ctxt total time : "<<ctxt_ctxt<<" s" << endl;
  cout<<"estimated TOTAL time : "<<ctxt_ctxt+hom_add+ptxt_ctxt<<" s" << endl;
  
  column=N_q*S;index = index%(N_q*S); // Delete this line for a real example
  vector<Plaintext> M_i(column); vector<vector<Plaintext>> M(N_a,M_i);

  for(int i=0;i<N_a;++i)
  {
    rand_diag_mat(M[i].begin(),M[i].size(),row, cc, lbd);
  }
  vector<int64_t> part(row);
  vector<vector<int64_t>> expected(N_a,part);
  for(int i=0;i<N_a;++i)
  {
    get_col_j(expected[i].begin(),M[i].begin(),index,column,row);
  }
  
  //vector<Ciphertext<DCRTPoly>> vect_i(N_q*S,E_pack);
  vector<vector<Ciphertext<DCRTPoly>>> all_extractions;//(cuts,vect_i);
  
  processingTime=0.0;
  TIC(t);
  
  full_extraction(all_extractions, E_Packs,lbd, S, cuts, cc);

  processingTime = TOC(t);
  
  cout << "extractions: " << processingTime/(1000) << "s"<<endl;
  
  vector<Ciphertext<DCRTPoly>> Answers(N_a,E_pack);
  
  full_evaluation(Answers, M.begin(), all_extractions.begin(), cuts, N_q*S, cc);
  
  for(int i=0;i<N_a;++i) 
  {
    if (!Serial::SerializeToFile(DATAFOLDER + "/answer_" +char(i)+ ".txt", cc->Compress(Answers[i]), SerType::BINARY)) 
    {
        std::cerr << "Error writting serialization of answer ctxt "<<i<< std::endl;
        return 1;
    }
    std::cout << "The answer ciphertexts have been serialized." << std::endl;
  }
  
  //VERIFICATION
  Plaintext testmat;
  
  for(int i=0;i<N_a;++i)
  {
  
    cc->Decrypt(kp.secretKey, Answers[i], &testmat);
    vector<int64_t> result = testmat->GetPackedValue();
    cout<<"VERIF "<<i<<" : "<<verification(expected[i], result,choice,row,cc)<<endl;
  }
  
  return 0;
}


