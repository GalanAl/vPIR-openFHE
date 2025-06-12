#include <iostream>
#include "openfhe.h"
#include <vector>
#include <algorithm>
#include <random>
#include "omp.h"

using namespace lbcrypto;
using namespace std;

const string DATAFOLDER = "SizeData";

//rotate inplace by k to the right
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


void fill_slots(Ciphertext<DCRTPoly>& V,const int32_t row, const CryptoContext<DCRTPoly>& cc)
{
  int32_t r=row;int count =0;long k=1;
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

// write index in corresponding base
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

void random_gen(vector<vector<int64_t>>& verif, const CryptoContext<DCRTPoly>& cc)
{
  srand (time(NULL));
  int64_t p = cc->GetCryptoParameters()->GetPlaintextModulus();
  for(int i=0; i<int(verif.size());++i) {for(long j=0;j<int(verif[i].size());++j) verif[i][j] = rand()%p;}  
  return;
}

void packing(Plaintext& pack, vector<bool>& choices, vector<vector<int64_t>>& verif, vector<vector<int64_t>>& query, const CryptoContext<DCRTPoly>& cc )
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  vector<int64_t> int_pack(slots);
  int64_t i;
  int lbd = choices.size(); 
  int64_t sep = slots/lbd;
  int L = query.size();
  int64_t S = query[0].size();
  for(int64_t k =0;k<lbd;++k)
  {
    i= k<lbd/2? k*sep: slots/2 + (k-lbd/2)*sep;
    if(choices[k])
    {
      for(int j=0;j<L;++j)
      {
        for(long f=0;f<S;++f) {int_pack[i]=query[j][f]; i+=1;}
      }
    }
    else
    {
      for(int j=0;j<L;++j)
      {
        for(long f=0;f<S;++f) {int_pack[i]=verif[j][f]; i+=1;}
      }
    }
  }
  
  pack = cc->MakePackedPlaintext(int_pack);
  
  processingTime = TOC(t);
  cout << "Packing time: " << processingTime << "ms" << std::endl;
  return;
}

void extract_one(vector<Ciphertext<DCRTPoly>>::iterator beg_V, const Ciphertext<DCRTPoly>& E_pack, const CryptoContext<DCRTPoly>& cc, const int32_t lbd, const int32_t S, const int num)
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  int32_t row = slots/lbd;
 
  for(long i=0;i<S;++i)
  {
    vector<int64_t> mask(slots);
    for (long l=0;l<lbd/2;++l) {mask[l*row+i+num*S]=1;mask[slots/2+ l*row+i+num*S]=1;}
    
    *(beg_V+i) = cc->EvalMult(E_pack , cc->MakePackedPlaintext(mask));
    
    rotate_k(*(beg_V+i),cc,-(i+S*num));
    
    fill_slots(*(beg_V+i),row,cc);
  }

  processingTime = TOC(t);
  cout << "Extract_one : " << processingTime/1000 << "s" << std::endl;
  return;
}

void extract_one_rot(vector<Ciphertext<DCRTPoly>>::iterator beg_V, const Ciphertext<DCRTPoly>& E_pack, const CryptoContext<DCRTPoly>& cc, const int32_t lbd, const int32_t S, const int num, const int L)
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  int32_t sep = slots/lbd;
  Ciphertext<DCRTPoly> temp;
  vector<int64_t> mask(slots),mask2(slots);
  for(long j=0;j<S;++j)
    {
      for (long l=0;l<lbd/2;++l) {mask[l*sep+j+num*S]=1;mask[slots/2+ l*sep+j+num*S]=1;}
    }
 
  for(long i=0;i<S;++i)
  {
    temp = cc->EvalMult(E_pack , cc->MakePackedPlaintext(mask2));
    rotate_k(temp,cc,-(S*num)+S-i);
    *(beg_V+i) = cc->EvalMult(E_pack , cc->MakePackedPlaintext(mask));
    rotate_k(*(beg_V+i),cc,-(S*num+i));
    *(beg_V+i) += temp;
    temp = *(beg_V+i);
    for(int k=1;k<L;++k){ rotate_k(temp,cc,S);*(beg_V+i) += temp;}
    
    for (long l=0;l<lbd/2;++l) {mask[l*sep+i+num*S]=0;mask[slots/2+ l*sep+i+num*S]=0;mask2[l*sep+i+num*S]=1;mask2[slots/2+ l*sep+i+num*S]=1;}
  }

  processingTime = TOC(t);
  cout << "Extract_one_rot : " << processingTime/1000 << "s" << std::endl;
  return;
}

void extract_coef(vector<Ciphertext<DCRTPoly>>::iterator beg_V, const Ciphertext<DCRTPoly>& E_pack, const CryptoContext<DCRTPoly>& cc, const int32_t lbd, const int64_t sep, const int32_t s_left, const int32_t s_right, const int32_t s_shift,const int32_t row)
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  #pragma omp parallel for shared(beg_V,E_pack,cc,lbd,s_left,s_right,sep,s_shift,slots,row)
  for(long i=0;i<s_left;++i)
  {
    vector<int64_t> mask(slots);
    for (long l=0;l<lbd/2;++l) {mask[l*sep+i+s_shift]=1;mask[slots/2+ l*sep+i+s_shift]=1;}
    
    Ciphertext<DCRTPoly> temp = cc->EvalMult(E_pack , cc->MakePackedPlaintext(mask));
    
    rotate_k(temp,cc,-(i+s_shift));
    Ciphertext<DCRTPoly> temp2;
    for(long k=1;k<row;k*=2)
    {
      temp2 = temp;
      rotate_k(temp2,cc,k);
      temp+= temp2;        
    }
    for(long j=0;j<s_right;++j)
    {
      *(beg_V+i*s_right+j) = temp;
    }
  }
  #pragma omp parallel for shared(beg_V,E_pack,cc,lbd,s_left,s_right,sep,s_shift,slots,row)
  for(long j=0;j<s_right;++j)
  {
    vector<int64_t> mask(slots);
    for (long l=0;l<lbd/2;++l) {mask[l*sep+j+s_shift+s_left]=1;mask[slots/2 +l*sep+j+s_shift+s_left]=1;}
    
    Ciphertext<DCRTPoly> temp = cc->EvalMult(E_pack , cc->MakePackedPlaintext(mask));
    
    rotate_k(temp,cc,-(j+s_shift+s_left));
    
    Ciphertext<DCRTPoly> temp2;
    for(long k=1;k<row;k*=2)
    {
      temp2 = temp;
      rotate_k(temp2,cc,k);
      temp+= temp2;        
    }
    for(long i=0;i<s_left;++i)
    {
      *(beg_V+i*s_right+j) *= temp;
    }
  }
  
  processingTime = TOC(t);
  cout << "Extract_coef : " << processingTime/1000 << "s" << std::endl;
  return;
}
void extract_rot(vector<Ciphertext<DCRTPoly>>::iterator beg_V, const Ciphertext<DCRTPoly>& E_pack, const CryptoContext<DCRTPoly>& cc, const int32_t lbd, const int64_t sep,const int32_t s_left, const int32_t s_right, const int32_t s_shift,const int32_t row)
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  vector<Ciphertext<DCRTPoly>> leftpart(s_left,*(beg_V));
  vector<Ciphertext<DCRTPoly>>::iterator beg_L = leftpart.begin();
  vector<Ciphertext<DCRTPoly>> rightpart(s_right,*(beg_V));
  vector<Ciphertext<DCRTPoly>>::iterator beg_R = rightpart.begin();
  
  #pragma omp parallel for shared(beg_L,E_pack,cc,lbd,sep,s_shift,slots,s_left,s_right)
  for(long i =0;i<s_left;++i)
  { 
    Ciphertext<DCRTPoly> temp;
    vector<int64_t> mask(slots);     
    
    for(long l =0;l<lbd/2;++l)
    {
      mask[l*sep+s_shift+i]=1;mask[slots/2+ l*sep+s_shift+i]=1;
    }

    *(beg_L+i) = cc->EvalMult(E_pack,cc->MakePackedPlaintext(mask));
      
    rotate_k(*(beg_L+i) ,cc,-(s_shift+i));
  
    for(long j=0;j<row/s_right;++j) // requires row+s_right < sep!!!  
    {
      mask.clear();mask.resize(slots);
      for(long l = 0;l<lbd/2;++l)
      {
        mask[l*sep+s_shift+((i+(j+1))%s_left)]=1;mask[slots/2+ l*sep+s_shift+((i+(j+1))%s_left)]=1;
      }
      
      temp=cc->EvalMult(E_pack,cc->MakePackedPlaintext(mask));
        
      rotate_k(temp,cc,-(s_shift+((i+(j+1))%s_left))+(j+1)*s_right);
      *(beg_L+i) += temp;
    }           
    for(long k=1;k<s_right;k*=2)
    {
      temp = *(beg_L+i);
      rotate_k(temp,cc,k);
      *(beg_L+i)+= temp;
    } 
  }    
  
  #pragma omp parallel for shared(beg_R,E_pack,cc,lbd,sep,s_shift,slots,s_left,s_right)
  for(long j=0;j<s_right;++j)
  {
    vector<int64_t> mask1(slots),mask2(slots); 
    Plaintext p_mask;
    Ciphertext<DCRTPoly> temp;
    
    for(long l=0;l<lbd/2;++l)
    {
      for(long i=0;i<j;++i)
      {
        mask1[s_shift+l*sep+s_left+i]=1;mask1[slots/2+s_shift+l*sep+s_left+i]=1;
      } 
      for(long i=j;i<s_right;++i)
      {
        mask2[s_shift+l*sep+s_left+i]=1;mask2[slots/2+s_shift+l*sep+s_left+i]=1;
      } 
    }

    *(beg_R+j) = cc->EvalMult(E_pack,cc->MakePackedPlaintext(mask2));
    rotate_k(*(beg_R+j),cc,-(s_shift+s_left+j));
    
    temp = cc->EvalMult(E_pack,cc->MakePackedPlaintext(mask1));
    rotate_k(temp,cc,-(s_shift+s_left)+s_right-j);
    
    *(beg_R+j) += temp;
    for(long k=1;k<row/s_right;k*=2) 
    {
      temp = *(beg_R+j);
      rotate_k(temp,cc,s_right*k);
      *(beg_R+j)+= temp;
    } 
  }
  
  #pragma omp parallel for shared(beg_V,beg_R,beg_L,s_left,s_right,cc)
  for(long j=0;j<s_left;++j)
  {
    for(long i = 0;i<s_right;++i)
    {
      *(beg_V+j*s_right+i) = cc->EvalMult(*(beg_R+i),*(beg_L+j));
      rotate_k(*(beg_L+j),cc,-1);
    }
  }
  processingTime = TOC(t);
  cout << "Extract_rot: " << processingTime/1000 << "s" << std::endl;
  return;
}

void rand_diag_mat(vector<vector<int64_t>>::iterator beg_M,const int64_t column, const CryptoContext<DCRTPoly>& cc, const int32_t lbd)
{
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  
  int64_t p = cc->GetCryptoParameters()->GetPlaintextModulus();
  int32_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
  int64_t row = slots/lbd;

  #pragma omp parallel for shared(beg_M,column,lbd,p,slots)
  for(long i = 0;i<column;++i)
  {
      srand (time(NULL));
      vector<int64_t> diag_i(slots);
      
      for(long j = 0;j<row;++j) 
      {
        int64_t r = rand()%p;
        for(long k = 0;k<lbd/2;++k) {diag_i[k*row+j]=r;diag_i[slots/2+k*row+j]=r;}
      }
      *(beg_M+i) = diag_i;
  }
  
  processingTime = TOC(t);
  cout << "matrix_gen: " << processingTime/1000 << "s" << std::endl;
  return;
}

void get_col_j(vector<int64_t>::iterator beg_col, vector<vector<int64_t>>::const_iterator beg_M, const int64_t j, const int64_t column, const int64_t row)
{
  vector<int64_t> temp;
  for(long i =0;i<row;++i)
  {
    temp = *(beg_M+((j-i)%column+column)%column);
    *(beg_col+i) = temp[i];
  }
  return;
}

void matrix_vector(Ciphertext<DCRTPoly>& MijV2, vector<vector<int64_t>>::const_iterator beg_M, vector<Ciphertext<DCRTPoly>>::const_iterator beg_rot, const CryptoContext<DCRTPoly>& cc, const int64_t column)
{
  /*
  int64_t sqr = sqrt(column);
  vector<Ciphertext<DCRTPoly>> baby(sqr),giant(sqr);
  for(long i=0;i<sqr;++i)
  {
    for(long j =0;j<sqr;++j)
    {
      baby[j] = cc->EvalMult(cc->MakePackedPlaintext(*(beg_M+sqr*i+j)),*(beg_rot+sqr*i+j));
    }
    giant[i] = cc->EvalAddMany(baby);
  }
  MijV2 = cc->EvalAddMany(giant);
  */
  vector<Ciphertext<DCRTPoly>> T(column);
  for(long i=0;i<column;++i) T[i] = cc->EvalMult(cc->MakePackedPlaintext(*(beg_M+i)),*(beg_rot+i));
  
  MijV2 = cc->EvalAddMany(T);
  return;
}

void full_evaluation(Ciphertext<DCRTPoly>& Answer, vector<vector<int64_t>>::const_iterator beg_M, vector<vector<Ciphertext<DCRTPoly>>>::const_iterator beg_Extr, const int L, const int32_t S, const CryptoContext<DCRTPoly>& cc)
{
  Ciphertext<DCRTPoly> T;
  switch(L)
  {
    case 1:
    {
      matrix_vector(Answer, beg_M, (*(beg_Extr)).begin(),cc,S); break;
    }
    case 2:
    {
      vector<Ciphertext<DCRTPoly>> A(S,T);
      vector<Ciphertext<DCRTPoly>>::iterator beg_A = A.begin();
      #pragma omp parallel for shared(cc,beg_A,beg_M,beg_Extr)
      for(int32_t a=0;a<S;++a)
      {
        matrix_vector(*(beg_A+a), beg_M, (*(beg_Extr+1)).begin(),cc,S);
        *(beg_A+a) = cc->EvalMult(*(beg_A+a),(*(beg_Extr))[a]);  
      }
      Answer = cc->EvalAddMany(A); break;
    }
    case 3:
    {
      int32_t size = pow(S,2);
      vector<Ciphertext<DCRTPoly>> AB(size,T);
      vector<Ciphertext<DCRTPoly>>::iterator beg_AB = AB.begin();
      #pragma omp parallel for shared(cc,beg_AB,S,beg_Extr,beg_M)
      for(int32_t i=0;i<size;++i)
      {
        matrix_vector(*(beg_AB+i), beg_M, (*(beg_Extr+2)).begin(),cc,S);
        *(beg_AB+i) = cc->EvalMult(*(beg_AB+i),cc->EvalMult((*(beg_Extr))[i/S],(*(beg_Extr+1))[i%S]));
      }
      cout<<"vect fill DONE"<<endl;
      Answer = cc->EvalAddMany(AB);break;   
    }
    case 4:
    {
      int32_t size = pow(S,2);
      vector<Ciphertext<DCRTPoly>> A(size,T);
      vector<Ciphertext<DCRTPoly>>::iterator beg_A = A.begin();
      #pragma omp parallel for shared(cc,beg_A,beg_Extr,beg_M,S)
      for(int32_t a=0;a<size;++a)
      {
        Ciphertext<DCRTPoly> Temp;
        vector<Ciphertext<DCRTPoly>> B(S,Temp);
        for(int32_t b=0;b<S;++b)
        {
          matrix_vector(Temp, beg_M, (*(beg_Extr+3)).begin(),cc,S);
          B[b] = cc->EvalMult(Temp,(*(beg_Extr+2))[b]);
        }
        *(beg_A+a)= cc->EvalMult(cc->EvalMult((*(beg_Extr))[a/S],(*(beg_Extr+1))[a%S]),cc->EvalAddMany(B));
      }
      Answer = cc->EvalAddMany(A);break;   
    }
    case 5:
    {
      int32_t size = pow(S,2);
      vector<Ciphertext<DCRTPoly>> A(size,T),B(size,T);
      vector<Ciphertext<DCRTPoly>>::iterator beg_A = A.begin();
      vector<Ciphertext<DCRTPoly>>::iterator beg_B = B.begin();
      #pragma omp parallel for shared(cc,beg_A,beg_B,beg_Extr,S)
      for(int32_t a=0;a<size;++a)
      {
        *(beg_B+a)=cc->EvalMult((*(beg_Extr+2))[a/S],(*(beg_Extr+3))[a%S]);
      }
      #pragma omp parallel for shared(cc,beg_A,beg_B,beg_Extr,S)
      for(int32_t a=0;a<size;++a)
      { 
        Ciphertext<DCRTPoly> Temp;
        vector<Ciphertext<DCRTPoly>> C(size,Temp);
        for(int32_t b=0;b<size;++b)
        {
          matrix_vector(Temp, beg_M, (*(beg_Extr+4)).begin(),cc,S);
          C[b]=cc->EvalMult(Temp,*(beg_B+b));
        }
        *(beg_A+a)=cc->EvalMult(cc->EvalAddMany(C),cc->EvalMult((*(beg_Extr))[a/S],(*(beg_Extr+1))[a%S]));
      }
      Answer = cc->EvalAddMany(A);break;   
    }
    case 6:
    {
      int32_t size = pow(S,2);
      vector<Ciphertext<DCRTPoly>> A(size,T),B(size,T);
      vector<Ciphertext<DCRTPoly>>::iterator beg_A = A.begin();
      vector<Ciphertext<DCRTPoly>>::iterator beg_B = B.begin();
      #pragma omp parallel for shared(cc,beg_A,beg_B,beg_Extr,S)
      for(int32_t a=0;a<size;++a)
      {
        *(beg_B+a)=cc->EvalMult((*(beg_Extr+2))[a/S],(*(beg_Extr+3))[a%S]);
      }
      #pragma omp parallel for shared(cc,beg_A,beg_B,beg_Extr,S)
      for(int32_t a=0;a<size;++a)
      { 
        Ciphertext<DCRTPoly> Temp;
        vector<Ciphertext<DCRTPoly>> C(size,Temp);
        vector<Ciphertext<DCRTPoly>> D(S,Temp);
        for(int32_t b=0;b<size;++b)
        {
          for(int32_t c=0;c<S;++c)
          {
            matrix_vector(Temp, beg_M, (*(beg_Extr+5)).begin(),cc,S);
            D[c]=cc->EvalMult((*(beg_Extr+4))[c],Temp);
          }
          C[b]=cc->EvalMult(cc->EvalAddMany(D),*(beg_B+b));
        }
        *(beg_A+a)=cc->EvalMult(cc->EvalAddMany(C),cc->EvalMult((*(beg_Extr))[a/S],(*(beg_Extr+1))[a%S]));
      }
      Answer = cc->EvalAddMany(A);break;  
    }
  }
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
        if(expected[j]!=result[(i>=lbd/2?slots/2+(i-lbd/2)*sep:i*sep)+j]){cout<<"PERDUUUU "<<i<<endl; return false;}
      }
    }
  }
  return true;
}

int main(int argc, char* argv[]) {  
    
    if(argc <= 1)
    {
      cerr << "Usage: " << argv[0] << " number of small pieces (max 6) " << endl;
      exit(0);
    }
    int L = (argc>1?atoi(argv[1]):6);
    
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

    parameters.SetPlaintextModulus(281474978414593);//
    parameters.SetMaxRelinSkDeg(3);
    parameters.SetMultiplicativeDepth(5);
    parameters.SetSecurityLevel(HEStd_128_classic);
    

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    // enable features that you wish to use
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    /*
    if (!Serial::SerializeToFile(DATAFOLDER + "/cryptocontext.txt", cc, SerType::BINARY)) {
        std::cerr << "Error writing serialization of the crypto context to "
                     "cryptocontext.txt"
                  << std::endl;
        return 1;
    }
    std::cout << "The cryptocontext has been serialized." << std::endl;
    */
    int64_t slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    cout<<"nb slot : "<<slots<<endl;
        
    // Initialize the public key containers.
    KeyPair<DCRTPoly> kp = cc->KeyGen();
    // 13 = log(slots)/2
    std::vector<int32_t> indexList; for(int32_t i=1;i<slots/2;i*=2){ indexList.push_back(int32_t(i));indexList.push_back(int32_t(-i));}
    cc->EvalRotateKeyGen(kp.secretKey, indexList);
    cc->EvalMultKeyGen(kp.secretKey);
    
    /*
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
  */
    //********************* CLIENT SIDE *********************//
  int32_t lbd = 42,sep = slots/lbd, S = sep/L, row = S*L;
  int64_t column = pow(S,L);
  cout<<"lbd = "<<lbd<<"; row = "<<row<<"; column = "<<column<<"; size vectors S = "<<S<<endl;
  
  vector<int64_t> length(S,0);
  vector<vector<int64_t>> query(L,length);
  vector<vector<int64_t>> verif(L,length);
  
  vector<bool> choice(lbd,true);
  srand (time(NULL));
  
  for(long i=0;i<lbd;++i) choice[i] = bool(rand()%2);
  choice[0]=1;
  
  int64_t index = rand()%int64_t(column);
  cout<<"index : "<<index<<" = "<<index%S<<" mod S"<<endl;
  query_gen(query, index);
  
  random_gen(verif, cc);
  
  /*
  cout<<"query : "<<query<<endl;
  cout<<"verif : "<<verif<<endl;
  cout<<"choice : "<<choice<<endl;
  */
  Plaintext pack;
  
  packing(pack, choice, verif, query, cc);
   
  Ciphertext<DCRTPoly> E_pack = cc->Encrypt(kp.publicKey, pack); 

  /*
  if (!Serial::SerializeToFile(DATAFOLDER + "/" + "firstmessage.txt", E_pack, SerType::BINARY)) {
        std::cerr << "Error writing serialization of firstmessage to firstmessage.txt" << std::endl;
        return 1;
    }
    std::cout << "The first ciphertext has been serialized." << std::endl;
  */
  vector<int64_t> diag(slots);
  
  // Change S for column for a real life example
  column=S;
  vector<vector<int64_t>> M(S,diag);vector<vector<int64_t>>::iterator beg_M = M.begin();
  
  rand_diag_mat(beg_M,M.size(), cc, lbd);
  vector<int64_t> expected(row);
  get_col_j(expected.begin(),beg_M,index%S,column,row);
  //cout<<"column "<<index%S<<" : "<<expected<<endl;
 
  Ciphertext<DCRTPoly> Answer;
  vector<Ciphertext<DCRTPoly>> E_V(S,Answer);
  vector<vector<Ciphertext<DCRTPoly>>> Extacted(L,E_V);
  
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  
  for(int i=0;i<L-1;++i) extract_one(Extacted[i].begin(), E_pack, cc, lbd, S, i);
  extract_one_rot(Extacted[L-1].begin(), E_pack,cc,lbd,S, L-1,L);
  
  processingTime = TOC(t);
  
  cout << "extractions: " << processingTime/(1000) << "s"<<endl;
  
  // TEST TIMINGS
  
  processingTime=0.0;
  TIC(t);
  
  #pragma omp parallel for shared(cc,Extacted,beg_M)
  for (long l=0;l<1000;++l)
  {
    Ciphertext<DCRTPoly> test_timings;
    matrix_vector(test_timings, beg_M, Extacted[L-1].begin(),cc,S);
    vector<Ciphertext<DCRTPoly>> Vec(L,test_timings);
    test_timings = cc->EvalMultMany(Vec);
    test_timings += test_timings;
  }
  processingTime = TOC(t);
  
  cout << "matrix vector: " << processingTime/(1000*1000) << "s; estimated time : "<<int64_t(pow(S,L-1)*processingTime*numthreads)/(56*1000*1000)<<" s" << std::endl;
  
  //FULL EVAL 
  processingTime=0.0;
  TIC(t);
  full_evaluation(Answer,beg_M, Extacted.begin(), L, S, cc);
  
  processingTime = TOC(t);
  
  cout << "process time: " << processingTime/(1000)<<"s"<<endl;
  
  Plaintext testmat;
  
  cc->Decrypt(kp.secretKey, Answer, &testmat);
  
  
  vector<int64_t> result = testmat->GetPackedValue();
  cout<<"VERIF: "<<verification(expected, result,choice,sep,cc)<<endl;

  //cout<<"result "<<result<<endl;
  
  return 0;
}

/*
  Ciphertext<DCRTPoly> none;
  vector<Ciphertext<DCRTPoly>> E_V0(s_0*s_1,none),E_V1(s_2*s_3,none),E_V2(s_4*s_5,none);
  vector<Ciphertext<DCRTPoly>>::iterator beg_V0= E_V0.begin(),beg_V1= E_V1.begin(),beg_V2= E_V2.begin();
  
  Plaintext testos;
  
  extract_coef(beg_V0, E_pack, cc, lbd, sep, s_0, s_1,0,row);
  
  extract_coef(beg_V1, E_pack, cc, lbd, sep, s_2, s_3, s_0+s_1,row);
  
  extract_rot(beg_V2, E_pack, cc, lbd, sep, s_4, s_5, s_0+s_1+s_2+s_3,row);
  
  VERIF MATVEC
  Ciphertext<DCRTPoly> E_testmat;
  Plaintext testmat;
  matrix_vector(E_testmat, beg_M, beg_V2,cc,column);
  
  cc->Decrypt(kp.secretKey, E_testmat, &testmat);
  cout<<"verif mat vec : "<<testmat<<endl;
  
  // TEST TIMINGS
  Ciphertext<DCRTPoly> test_timings;
  TimeVar t;
  double processingTime(0.0);
  TIC(t);
  for (long l=0;l<100;++l)
  {
    matrix_vector(test_timings, beg_M, beg_V2,cc,column);
  }
  processingTime = TOC(t);
  
  cout << "matrix vector: " << processingTime/(1000*100) << "s" << std::endl;
  
  processingTime=0.0;
  TIC(t);
  for (long l=0;l<100;++l)
  {
    test_timings = cc->EvalMult(*(beg_V1),*(beg_V0));
  }
  processingTime = TOC(t);
  
  cout << "first mult: " << processingTime/(1000*100) << "s" << std::endl;
  
  Ciphertext<DCRTPoly> test_timings2;
  processingTime=0.0;
  TIC(t);
  for (long l=0;l<100;++l)
  {
    test_timings2 = cc->EvalMult(test_timings,test_timings);
  }
  processingTime = TOC(t);
  
  cout << "second mult: " << processingTime/(1000*100) << "s" << std::endl;
  
  Ciphertext<DCRTPoly> full_sum;
  processingTime=0.0;
  TIC(t);
  full_sum = test_timings2;
  for (long l=1;l<pow(2,20);++l)
  {
    full_sum += test_timings2;
  }
  processingTime = TOC(t);
  
  cout << "2**20 additions: " << processingTime/(1000) << "s" << std::endl;
  
  
  Ciphertext<DCRTPoly> E_M;
  vector<Ciphertext<DCRTPoly>> Temp(numthreads,E_M);vector<Ciphertext<DCRTPoly>>::iterator beg_Temp = Temp.begin();
  vector<Ciphertext<DCRTPoly>> Temp2(numthreads,E_M);vector<Ciphertext<DCRTPoly>>::iterator beg_Temp2 = Temp2.begin();
  int32_t S_0=s_0*s_1;  
  
  processingTime=0.0;
  TIC(t);
  
  //#pragma omp parallel for shared(cc,beg_Temp,beg_V0,beg_V1,beg_V2,column,beg_M)
  for(long i=0;i<numthreads;++i)
  {
    matrix_vector(*(beg_Temp+i), beg_M, beg_V2,cc,column);
    *(beg_Temp+i) *= cc->EvalMult(*(beg_V1),*(beg_V0+i));
  }
  processingTime = TOC(t);
  
  cout << "Thing to do many times: " << processingTime/(1000*numthreads) << "s" << std::endl;  
    
  processingTime=0.0;
  TIC(t);  
    
  #pragma omp parallel for shared(cc,beg_Temp,beg_V0,beg_V1,beg_V2,column,beg_M,S_0,beg_Temp2)
  for(long i=numthreads;i<S_0;++i)
  {
    matrix_vector(*(beg_Temp2+omp_get_thread_num()), beg_M, beg_V2,cc,column);
    *(beg_Temp2+omp_get_thread_num()) *= cc->EvalMult(*(beg_V1),*(beg_V0+i));
    *(beg_Temp+omp_get_thread_num()) += *(beg_Temp2+omp_get_thread_num());
  }
  cout<<"THE END"<<endl;
  
  E_M = cc->EvalAddManyInPlace(Temp);
  
  processingTime = TOC(t);
  
  cout << "partial eval: " << processingTime/1000 << "s" << std::endl;
  
  if (!Serial::SerializeToFile(DATAFOLDER + "/" + "secondmessage.txt", E_M, SerType::BINARY)) {
        std::cerr << "Error writing serialization of secondmessage to secondmessage.txt" << std::endl;
        return 1;
    }
    std::cout << "The second ciphertext has been serialized." << std::endl;

  cc->Decrypt(kp.secretKey, E_M, &testos);
  cout<<"final verif : "<<testos<<endl;*/

