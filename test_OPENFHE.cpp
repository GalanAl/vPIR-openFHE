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

void full_extraction(vector<vector<Ciphertext<DCRTPoly>>>::iterator beg_extract, const vector<Ciphertext<DCRTPoly>>& E_Packs, const int lbd, const int32_t S, const int cuts, const CryptoContext<DCRTPoly>& cc)
{
  for(int i = 0;i<cuts-1;++i)
  {
    extract_one((*(beg_extract+i)).begin(),E_Packs,lbd, S, i, S*cuts, cc);
  }
  extract_one_rot((*(beg_extract+cuts-1)).begin(), E_Packs,lbd, S,cuts,cc);
  return;
}

void verif_extraction(vector<vector<Ciphertext<DCRTPoly>>>& All_extract, vector<bool>& choices, vector<vector<int64_t>>& verif, vector<vector<int64_t>>& query,const int32_t row, const CryptoContext<DCRTPoly>& cc,KeyPair<DCRTPoly>& kp)
{
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  Plaintext Ptest;
  vector<int64_t> test;
  int32_t jump;
  int64_t value_q,value_v;
  int cuts = int(All_extract.size());
  for(int i=0;i<cuts-1;++i)
  {
    for(int32_t j=0;j<int32_t(All_extract[i].size());++j)
    {
      cc->Decrypt(kp.secretKey, All_extract[i][j], &Ptest);
      test = Ptest->GetPackedValue();
      value_q=query[i][j];
      value_v=verif[i][j];
      for(int k=0;k<int(choices.size());++k)
      {
        jump = k<int(choices.size()+1)/2?0:slots/2-row*((choices.size()+1)/2);
        for(int32_t l=0;l<row;++l)
        {
          if(test[jump+k*row+l]!=(choices[k]?value_q:value_v))
          {
            cout<<"FAILURE i="<<i<<" ; j="<<j<<" ; k="<<k<<" ; l="<<l<<endl;
            return;
          }
        }
      }
    } 
  }
  for(int32_t j=0;j<int32_t(All_extract[cuts-1].size());++j)
  {
    cc->Decrypt(kp.secretKey, All_extract[cuts-1][j], &Ptest);
    test = Ptest->GetPackedValue();
    for(int k=0;k<int(choices.size());++k)
    {
      jump = k<int(choices.size()+1)/2?0:slots/2-row*((choices.size()+1)/2);
      for(int32_t l=0;l<row;++l)
      {
        if(test[jump+k*row+l]!=(choices[k]?query[cuts-1][(j+l)%query[cuts-1].size()]:verif[cuts-1][(j+l)%verif[cuts-1].size()]))
        {
          cout<<"FAILURE i="<<cuts-1<<" ; j="<<j<<" ; k="<<k<<" ; l="<<l<<endl;
          return;
        }
      }
    }
  }
}

// generate a random matrix (arranged by diagonals)
void rand_diag_mat(vector<vector<int64_t>>::iterator beg_M,const int64_t column,const int32_t row, const CryptoContext<DCRTPoly>& cc, const int lbd)
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
      *(beg_M+i) = diag_i;
  }
  
  processingTime = TOC(t);
  cout << "matrix_gen: " << processingTime/1000 << "s" << std::endl;
  return;
}

// get the j^th column of a matrix described by its diagonals
void get_col_j(vector<int64_t>::iterator beg_col, vector<vector<int64_t>>::const_iterator beg_M, const int64_t j, const int64_t column, const int64_t row)
{
  vector<int64_t> temp;
  for(long i =0;i<row;++i)
  {
    temp = *(beg_M+((j-i)%column+column)%column);//TODO??!!!
    *(beg_col+i) = temp[i];
  }
  return;
}

// matrix-vector product ; sum of each diagonals times rotations of the vector
void matrix_vector(Ciphertext<DCRTPoly>& Mat_vec, vector<vector<int64_t>>::const_iterator beg_M_i, vector<Ciphertext<DCRTPoly>>::const_iterator beg_rot, const CryptoContext<DCRTPoly>& cc, const int32_t column)
{/*
  int32_t sqr = sqrt(column)+1;
  vector<Ciphertext<DCRTPoly>> giant(sqr);
  int32_t j= 0;
  for(int32_t i=0;i<sqr;++i)
  {
    vector<Ciphertext<DCRTPoly>> baby(min(column-i*sqr,sqr));
    for(;j<min(column,(i+1)*sqr);++j)
    {
      baby[j%sqr] = cc->EvalMult(cc->MakePackedPlaintext(*(beg_M_i+j)),*(beg_rot+j));
    }
    giant[i] = cc->EvalAddMany(baby);
  }
  Mat_vec = cc->EvalAddMany(giant);
  */
  vector<Ciphertext<DCRTPoly>> T(column,Mat_vec);  
  for(int32_t i=0;i<column;++i)
  {
    T[i] = cc->EvalMult(cc->MakePackedPlaintext(*(beg_M_i+i)),*(beg_rot+i));
  }
  Mat_vec = cc->EvalAddMany(T);
  
  return;
}

void full_evaluation(vector<Ciphertext<DCRTPoly>>& Answers, vector<vector<vector<int64_t>>>::const_iterator beg_M, vector<vector<Ciphertext<DCRTPoly>>>::const_iterator beg_extract, const int cuts, const int32_t column, const CryptoContext<DCRTPoly>& cc)
{
  int64_t numthreads(1);
  #if defined(_OPENMP)
    #pragma omp parallel
    #pragma omp single
    {
      numthreads = omp_get_num_threads();
    }
    #endif
    
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  

  switch(cuts)
  {
    case 1:
    {
      for(int i=0;i<int(Answers.size());++i) matrix_vector(Answers[i], (*(beg_M+i)).begin(), (*(beg_extract)).begin(),cc,column);
      break;
    }
    case 2:
    {
      for(int i=0;i<int(Answers.size());++i) 
      {
        vector<Ciphertext<DCRTPoly>> Accu(column,Answers[i]);
        vector<Ciphertext<DCRTPoly>>::iterator beg_A = Accu.begin();
        #pragma omp parallel for shared(beg_A,beg_M,beg_extract)
        for(int32_t A=0;A<column;++A)
        {
          matrix_vector(*(beg_A+A), (*(beg_M+i)).begin(), (*(beg_extract+1)).begin(),cc,column);//+A for second argument in real life 
          (*(beg_A+A)) *= (*(beg_extract))[A];
        }
        Answers[i] = cc->EvalAddMany(Accu);
      }
      break;
    }
    case 3:
    {
      for(int i=0;i<int(Answers.size());++i) 
      {
        vector<Ciphertext<DCRTPoly>> Accu(column,Answers[i]);
        vector<Ciphertext<DCRTPoly>>::iterator beg_A = Accu.begin();
        
        #pragma omp parallel for shared(beg_A,beg_M,beg_extract)
        for(int32_t A=0;A<column;++A)
        { 
          vector<Ciphertext<DCRTPoly>> Bccu(column,Answers[i]);
          for(int32_t B=0;B<column;++B)
          {
            matrix_vector(Bccu[B], (*(beg_M+i)).begin(), (*(beg_extract+2)).begin(),cc,column);//+A*column+B for second argument in real life 
            Bccu[B] *= cc->EvalMult((*(beg_extract))[A],(*(beg_extract+1))[B]);
          }
          (*(beg_A+A)) = cc->EvalAddMany(Bccu);
        }
        Answers[i] = cc->EvalAddMany(Accu);
      }
      break;
    }
    case 4:
    {
      for(int i=0;i<int(Answers.size());++i) 
      {
        vector<Ciphertext<DCRTPoly>> Accu(column,Answers[i]);
        vector<Ciphertext<DCRTPoly>>::iterator beg_A = Accu.begin();
        #pragma omp parallel for shared(beg_A,beg_M,beg_extract)
        for(int32_t A=0;A<column;++A)
        {
          vector<Ciphertext<DCRTPoly>> Bccu(column,*beg_A);
          vector<Ciphertext<DCRTPoly>> Cccu(column,*beg_A);
          for(int32_t B=0;B<column;++B)
          {
            for(int32_t C=0;C<column;++C)
            {
              matrix_vector(Cccu[C], (*(beg_M+i)).begin(), (*(beg_extract+3)).begin(),cc,column);//+A*column**2+B*column+C for second argument in real life 
              Cccu[C] *= (*(beg_extract+2))[C];
              Cccu[C] *= cc->EvalMult((*(beg_extract))[A],(*(beg_extract+1))[B]);
            }
            Bccu[B] = cc->EvalAddMany(Cccu);
          }
          (*(beg_A+A)) = cc->EvalAddMany(Bccu);
        }
        Answers[i] = cc->EvalAddMany(Accu);
      }
      break;
    }
    case 5:
    {
      for(int i=0;i<int(Answers.size());++i) 
      {
        vector<Ciphertext<DCRTPoly>> Accu(column,Answers[i]);
        vector<Ciphertext<DCRTPoly>>::iterator beg_A = Accu.begin();
        #pragma omp parallel for shared(beg_A,beg_M,beg_extract)
        for(int32_t A=0;A<column;++A)
        {
          vector<Ciphertext<DCRTPoly>> Bccu(column,*beg_A);
          vector<Ciphertext<DCRTPoly>> Cccu(column,*beg_A);
          vector<Ciphertext<DCRTPoly>> Dccu(column,*beg_A);
          for(int32_t B=0;B<column;++B)
          {
            for(int32_t C=0;C<column;++C)
            {
              for(int32_t D=0;D<column;++D)
              {
                matrix_vector(Dccu[D], (*(beg_M+i)).begin(), (*(beg_extract+4)).begin(),cc,column);//+A*column**3+B*column**2+C*column+D for second argument in real life 
                Dccu[D] *= cc->EvalMult(cc->EvalMult((*(beg_extract))[A],(*(beg_extract+1))[B]),cc->EvalMult((*(beg_extract+2))[C],(*(beg_extract+3))[D]));
              }
              Cccu[C] = cc->EvalAddMany(Dccu);
            }
            Bccu[B] = cc->EvalAddMany(Cccu);
          }
          (*(beg_A+A)) = cc->EvalAddMany(Bccu);
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
  
  Plaintext pack;
  vector<Plaintext> Packs(N_q,pack);
  Ciphertext<DCRTPoly> E_pack;
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
  
  
  column=N_q*S;index = index%(N_q*S); // Delete this line for a real example
  vector<vector<int64_t>> M_i(column,diag); vector<vector<vector<int64_t>>> M(N_a,M_i);

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
  
  //cout<<"EXPECTED "<<expected<<endl;
  
  vector<Ciphertext<DCRTPoly>> vect_i(N_q*S,E_pack);
  vector<vector<Ciphertext<DCRTPoly>>> all_extractions(cuts,vect_i);
  
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  
  full_extraction(all_extractions.begin(), E_Packs,lbd, S, cuts, cc);

  processingTime = TOC(t);
  
  cout << "extractions: " << processingTime/(1000) << "s"<<endl;
  
  #pragma omp parallel for
  for (long l=0;l<1000;++l)
  {
    Ciphertext<DCRTPoly> time;
    matrix_vector(time, M[0].begin(), all_extractions[cuts-1].begin(),cc,column);
    vector<Ciphertext<DCRTPoly>> multi(cuts-1,all_extractions[0][l%column]);
    time *= cc->EvalMultMany(multi);
  }
  processingTime = TOC(t);
  
  cout << "average time per element : " << processingTime/1000 << "ms; estimated total time: "<<(pow(column,cuts-1)*(processingTime/1000)*N_a)/1000<<" s" << std::endl;
  
  vector<Ciphertext<DCRTPoly>> Answers(N_a,E_pack);
  
  full_evaluation(Answers, M.begin(), all_extractions.begin(),cuts, N_q*S, cc);
  
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
    //cout<<"Expected ::: "<<expected[i]<<endl;
    //cout<<"Result ::: "<<result<<endl;
  }
  
  return 0;
}


