//
//   TTR - Tadpoles Tensor Reduction
//
//   Generator of FORM tensor reduction tables for tadpole integrals
//
//   Based on ideas from [arXiv:1701.014], and use FireFly for matrix inversion
//
//   Andrey Pikelner, 2023
//
//
//   Usage: echo "p1(mu1)*p2(mu2)..." | ./TTR [options]
// 
// 
#include "firefly/DenseSolver.hpp"
#include "firefly/Reconstructor.hpp"
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <deque>
#include <unordered_map>
#include <sstream>


using std::cout;
using std::cin;
using std::string;
using std::ostringstream;
using std::endl;
using std::make_pair;
using std::map;
using std::multimap;
using std::pair;
using std::set;
using std::vector;
using std::unordered_map;

// Polynomial with integer coefficients
// result of index contractions
typedef vector<unsigned long> DPoly;
typedef uint8_t SmallIdx;

// list of indices paired into metric tensors arguments
typedef std::pair<SmallIdx, SmallIdx> IdxPair;
typedef std::list<IdxPair> TensPairs;
// typedef std::deque<IdxPair> TensPairs;
typedef vector<size_t> LoopInfo;
typedef vector<SmallIdx> Perm;

// Finite field matrix inversion
namespace firefly
{
  template<typename FFIntTemp>
  FFIntTemp poly2ff (const DPoly & in_poly, const FFIntTemp& d)
  {
    FFIntTemp acc;
    for (int i = 0; i < in_poly.size(); ++i)
      {
        auto cn_ff = in_poly[i];
        if(cn_ff != 0 )
          {
            FFIntTemp c_ff(cn_ff);
            FFIntTemp pp(i);
            acc += c_ff*d.pow(pp);
          }
      }

    return acc;
  }


  template<typename FFIntTemp>
  void mat2ff(const vector<DPoly> & m_int, size_t n, mat_ff<FFIntTemp> & mff, const FFIntTemp& d)
  {
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
        mff[i][j] = poly2ff(m_int[n*i + j], d);
  }

  class BlackBoxUser : public BlackBoxBase<BlackBoxUser>
  {
  private:
    vector<DPoly> msymb;
    // Matrix size
    size_t n;
  public:
    BlackBoxUser(const vector<DPoly> &msymb_, size_t n_) : msymb(msymb_), n(n_) { };
    
    template<typename FFIntTemp>
    std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values)
    {
      
      FFIntTemp d(values.front());
      
      // Get results from parsed expressions
      std::vector<FFIntTemp> result;

      // Reserve (n x n) matrix 
      std::vector<std::vector<FFIntTemp>> mat(n,std::vector<FFIntTemp>(n));

      mat2ff(msymb, n, mat, d);
      calc_inverse(mat,mat.size());
      for (auto mrow : mat)
        for (auto m_ij : mrow)
          result.emplace_back(m_ij);

      return result;
    }

    inline void prime_changed() { }

  };

}

vector<firefly::RationalFunction> ff_inverse_and_select(const vector<DPoly> &m, size_t n, size_t col_pos, size_t wrkrs, bool quiet)
{
  auto verb = firefly::Reconstructor<firefly::BlackBoxUser>::IMPORTANT;
  if (quiet)
    verb = firefly::Reconstructor<firefly::BlackBoxUser>::SILENT;
  firefly::BlackBoxUser m_sampler(m, n);
  firefly::Reconstructor<firefly::BlackBoxUser>
    m_reconstructed(1, wrkrs,
                    m_sampler,
                    verb
                    // firefly::Reconstructor<firefly::BlackBoxUser>::SILENT
                    // firefly::Reconstructor<firefly::BlackBoxUser>::IMPORTANT
                    // firefly::Reconstructor<firefly::BlackBoxUser>::CHATTY
                    );
  m_reconstructed.reconstruct();
  std::vector<firefly::RationalFunction> results = m_reconstructed.get_result();

  std::vector<firefly::RationalFunction> sel_res;

  for (size_t i = 0; i < n; i++)
    sel_res.emplace_back(results[i*n + col_pos]);

  return sel_res;
}


// Version with non-paired indices in each tensor
LoopInfo contract(const Perm &t1, const Perm &t2)
{
  // TensPairs ct;
  TensPairs ct(t1.size());
  size_t i = 0;
  // Important, fill it sorted
  SmallIdx mu,nu;
  TensPairs::iterator it = ct.begin();
  for (Perm::const_iterator pit = t1.begin(); pit != t1.end();)
    {
      mu = *(pit++);
      nu = *(pit++);
      
      if (mu < nu)
        {
          it->first = mu;
          it->second = nu;
        }
      else
        {
          it->first = nu;
          it->second = mu;          
        }
      it++;
    }
  for (Perm::const_iterator pit = t2.begin(); pit != t2.end();)
    {
      mu = *(pit++);
      nu = *(pit++);
      
      if (mu < nu)
        {
          it->first = mu;
          it->second = nu;
        }
      else
        {
          it->first = nu;
          it->second = mu;          
        }
      it++;
    }

    // sort pairs by first index
  ct.sort();


  // Actual contraction
  size_t d_pow = 0;
  LoopInfo loops;
  size_t loop_len = 1;
  while (ct.size() > 0)
    {
      // We pop first elelment of the list, than modify elelemnt it can be contracted with
      // If we pop g_xx we simply remove it and increase n in d^n
      std::tie(mu, nu) = ct.front();
      // cout << "mu = " << mu << ", nu = " << nu << endl;
      if (mu == nu)
        {
          ct.pop_front();
          loops.push_back(loop_len / 2);
          loop_len = 1;
        }
      else
        {
          // contract nu, search starts from second element
          auto git = std::find_if(std::next(ct.begin()), ct.end(), [&](const IdxPair ab)
          {
            return (nu == ab.first) || (nu == ab.second);
          });
          
          // TODO
          if (git == ct.end())
            throw std::logic_error("Not found index during tensor contraction");
          
          
          ct.begin()->second = (git->first == nu) ? git->second : git->first;
          
          ct.erase(git);
          
          loop_len++;
        }
    }
  
  // Make representation unique
  std::sort(loops.begin(), loops.end());
  return loops;
}

// Generate vector 0...rank-1
Perm gen_zero_perm(size_t rank)
{
  Perm ivec(rank);
  std::iota(ivec.begin(), ivec.end(), 0);
  std::sort(ivec.begin(), ivec.end());
  return ivec;
}


std::ostream &operator<<(std::ostream &os, Perm const &p) {
  os << "<";
  for (auto pit = p.begin(); pit != p.end(); ++pit)
    if (pit == p.begin())
      os << " " << static_cast<unsigned int>(*pit);
    else
      os << ", " << static_cast<unsigned int>(*pit);
  os << " >";
  return os;
}

std::ostream &operator<<(std::ostream &os, LoopInfo const &p) {
  os << "{ ";
  for (auto pit = p.begin(); pit != p.end(); ++pit)
    os << static_cast<unsigned int>(*pit);
  os << " }";
  return os;
}


// All manipulations with entered numerator
class NumeratorTensor
{
  vector<string> vecs;
  vector<string> idxs;

public:
  NumeratorTensor(const std::string & input)
  {
    // parse product of vectors
    std::stringstream ss (input);
    std::string item;
    
    while (getline (ss, item,'*')) {
      item.erase(std::remove_if(item.begin(), item.end(), ::isspace), item.end());

      // k1(mu1)
      // position of "("
      std::size_t pos_obr = item.find("(");
      if (pos_obr == std::string::npos)
        throw std::logic_error("wrong momentum input");
      
      // create vector
      vecs.push_back(item.substr(0, pos_obr));
      
      // position of ")"
      std::size_t pos_cbr = item.find(")",pos_obr);
      if (pos_cbr == std::string::npos)
        throw std::logic_error("wrong indices input");        
      
      // create index
      idxs.push_back(item.substr(pos_obr+1,pos_cbr-pos_obr-1));

    }
    
  }

  size_t rank() const
  {
    return idxs.size();
  }

  // input info
  string info() const
  {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < vecs.size(); i++)
      {
        if (i > 0) ss << " ";
        ss << vecs[i] << "." << idxs[i];
      }
    ss << "]";
    return ss.str();
  }

  // LHS of id statement
  string id_lhs() const
  {
    std::stringstream ss;
    for (size_t i = 0; i < vecs.size(); i++)
      {
        if (i > 0) ss << "*";
        ss << vecs[i] << "(" << idxs[i] << "?)";
      }
    return ss.str();
  }

  // string of scalar product after contraction with metric tensors permutation
  string contract(const Perm & Tperm) const
  {
    std::stringstream ss;

    vector<string> v_sp;
    for (size_t i2 = 0; i2 < Tperm.size(); i2 += 2)
      {
        if (vecs[Tperm[i2]] < vecs[Tperm[i2+1]])
          v_sp.push_back(vecs[Tperm[i2]] + "." + vecs[Tperm[i2+1]]);
        else
          v_sp.push_back(vecs[Tperm[i2+1]] + "." + vecs[Tperm[i2]]);
      }
    std::sort(v_sp.begin(), v_sp.end());

    // Output unique representation
    for (auto sp_it = v_sp.begin(); sp_it != v_sp.end(); ++ sp_it)
      {
        if (sp_it != v_sp.begin()) ss << "*";
        ss << *sp_it;
      }
    return ss.str();
  }

  // string of metric tensors with original tensor indices permuted with Tperm
  string g_perm(const Perm & Tperm) const
  {
    std::stringstream ss;
    for (size_t i2 = 0; i2 < Tperm.size(); i2 += 2)
      {
        if(i2 > 0) ss << "*";
        ss << "d_(" << idxs[Tperm[i2]] << "," << idxs[Tperm[i2+1]] << ")";
      }
    return ss.str();
    
  }
};


class Permutations
{
  size_t rank;
  size_t max_loops;             // rank/2
  std::vector<std::vector<SmallIdx>> perms;

  // Minimal set of permutations giving different contraction
  // results after contraction with zero permutation
  map<LoopInfo, Perm> min_perms;

  // Multiplicity of each permutation
  vector<size_t> min_multiplicity;

  // Results of contraction between minimal permutations
  // Polynomial stored as [c0,c1,c2,c3,...] = c0 + d*c1 + d^2*c2 + ...
  // LoopInfo encodes result of contraction with p0
  map<pair<LoopInfo, LoopInfo>, DPoly> min_contr_polys;

  // Zero permutation to clasify all permutation
  Perm p0;

  // Size of minimal system to invert
  size_t m;

  // Solution for coefficients after matrix inversion
  vector<firefly::RationalFunction> b;

  // flags
  size_t workers;
  bool be_quiet;

public:
  using I = uint32_t;
  static constexpr size_t MAXRANK = sizeof(I) * __CHAR_BIT__;
  typedef std::bitset<MAXRANK> IntBits;

  Permutations(size_t rank_, size_t workers_ = 1, bool be_quiet_ = false);

  // Sort permutations according to the
  // contraction loop structure
  void loops();
  void solveB();

  // Export FORM file
  void form_export(std::ostream &, const NumeratorTensor&, bool);
};


Permutations::I bit_twiddle_permute(Permutations::I v)
{
  Permutations::I t = v | (v - 1); // t gets v's least significant 0 bits set to 1
  // Next set to 1 the most significant bit to change,
  // set to 0 the least significant ones, and add the necessary 1 bits.
  Permutations::I w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));
  return w;
}

// split vector v=v1,v2,v3,.... according to the bit mask
// into vA and vB with vA indices on 1s positions and rest sorted
template <typename T>
void split2(const std::vector<T> &v, const Permutations::IntBits &patt,
            std::vector<T> &v1, std::set<T> &v0)
{
  size_t pos1 = 0;
  size_t pos0 = 0;
  for (size_t pos = 0; pos < v.size(); pos++)
    if (patt.test(pos)) v1[pos1++] = v[pos];
    else v0.insert(v[pos]);
}


// glue two sets n1,n2,n3,n4,... + m1,m2,m3,m4,... -> (n1,m1)(n2,m2)...
template <typename T>
std::vector<T> createT(const std::vector<T> &v1, const std::vector<T> &v0)
{
  std::vector<T> v(2 * v1.size());
  for (size_t i = 0; i < v1.size(); i++)
    {
      v[2 * i] = (v1[i]);
      v[2 * i + 1] = (v0[i]);
    }
  return v;
}



// recursive function creating possible lists for sorted permutations
template <typename T>
void fill_nest(vector<vector<T>> &acc, const vector<T> &etalon,
               const vector<T> &sorted, const set<T> &rest)
{
  // Nothing to do, returning checked set
  if (etalon.size() == sorted.size())
    acc.push_back(sorted);
  else
    {
      // all rest elements must be greater than newely analysed "etalon" element,
      // since they are sorted and first indices are also sorted
      auto it = rest.upper_bound(etalon[sorted.size()]);
      if (it == rest.begin())
        for (; it != rest.end(); ++it)
          {
            // TODO improve
            auto lv = sorted;
            lv.push_back(*it);
            auto rs = rest;
            rs.erase(*it);
            
            // recursive call
            fill_nest(acc, etalon, lv, rs);
          }
    }
}

Permutations::Permutations(size_t rank_, size_t workers_, bool be_quiet_) : rank(rank_), workers(workers_), be_quiet(be_quiet_)
{
  Perm ivec(rank);
  std::iota(ivec.begin(), ivec.end(), 0);
  std::sort(ivec.begin(), ivec.end());

  // 00000...r/2|11111...r/2
  I p = pow(2, rank / 2) - 1;
  I n = p;
  do {
    Perm vA(rank / 2);
    std::set<SmallIdx> vB; // (rank/2);
    split2(ivec, IntBits(n), vA, vB);

    Perm srt;
    vector<Perm> acc;
    fill_nest(acc, vA, srt, vB);
      
    for (Perm aa : acc)
      perms.push_back(createT(vA, aa));
    // looping
    p = n;
    n = bit_twiddle_permute(p);
  } while (n > p && n < pow(2, rank) - 1);

  // All polynomials have maximal degree
  max_loops = rank / 2;

  // 
  loops();
  // Find solution for b_i
  solveB();

}


void Permutations::loops()
{
  // fix etalon perm
  p0 = gen_zero_perm(rank);
  
  // Calculate contractions with p0
  multimap<LoopInfo, Perm> sperm;
  for (auto it1 = perms.cbegin(); it1 != perms.cend(); ++it1)
    sperm.insert(make_pair(contract(p0, *it1), *it1));
  
  // Representative permutations for each loop configuration
  for (auto repr_it = sperm.begin(); repr_it != sperm.end();
       repr_it = sperm.upper_bound(repr_it->first))
    min_perms[repr_it->first] = repr_it->second;

  m = min_perms.size();
  
  cout << "Perms size (n) = " << perms.size() << endl
       << "Rank       (r) = " << rank << endl
       << "Vars num   (m) = " << m << endl;

  // iterate over representative permutations to contract with
  for (auto contr_it = min_perms.begin(); contr_it != min_perms.end(); ++contr_it)
    {
      // Now iterate over first elements of representative sets
      for (auto repr_it = sperm.begin(); repr_it != sperm.end();
           repr_it = sperm.upper_bound(repr_it->first))
        {          
          DPoly dp(rank/2 + 1, 0);
          // Now contract with all elements of the set
          for (auto eq_it = repr_it; eq_it != sperm.upper_bound(repr_it->first); ++eq_it)
            dp[contract(contr_it->second, eq_it->second).size()]++;

          min_contr_polys[make_pair(repr_it->first, contr_it->first)] = dp;
        }
      // Store muiltiplicity info
      min_multiplicity.push_back(sperm.count(contr_it->first));
    }
}


void Permutations::solveB()
{
  // Represent martix as vector
  vector<DPoly> mat_d(m*m);
  // Store position of 1 in rhs, to extract from unverse matrix
  set<size_t> rhs_one_pos;

  size_t i = 0;
  for (auto pit1 = min_perms.begin(); pit1 != min_perms.end(); ++pit1, i++)
    {
      size_t j = 0;
      for (auto pit2 = min_perms.begin(); pit2 != min_perms.end(); ++pit2, j++)
        {
          auto eit = min_contr_polys.find(make_pair(pit1->first, pit2->first));
          if (eit != min_contr_polys.end())
            mat_d[m*i + j] = eit->second;
          else
            throw std::logic_error("Not found coefficient");
        }
      
      // Assumed that we have single 1 in first position
      if (contract(p0, pit1->second).size() == max_loops) rhs_one_pos.insert(i);
      
    }
  
  // TODO
  if(rhs_one_pos.size() != 1)
    throw std::logic_error("More than one \"1\" in RHS");
  
  // Inverse matrix an return only needed column
  b = ff_inverse_and_select(mat_d, m, *(rhs_one_pos.begin()), workers, be_quiet);

  cout << " - Matrix to fix B_i inverted" << endl;
}


// Convert num,den to polyratfun
string nd2rat(const firefly::RationalFunction & r)
{
  std::vector<firefly::Monomial> num = r.numerator.coefs;
  std::vector<firefly::Monomial> den = r.denominator.coefs;

  std::stringstream ss;
  ss << "rat(";
  // numerator
  for (vector<firefly::Monomial>::const_iterator n_it = num.begin(); n_it != num.end(); ++n_it)
    {
      
      // TODO check zero
      ss << "+ ("
         << n_it->coef.numerator.get_str() <<  "/"
         << n_it->coef.denominator.get_str() << ")";
      if (n_it->powers.front() > 0)
        ss << "*d^" << n_it->powers.front();
    }
  ss << ",";
  // denominator
  for (vector<firefly::Monomial>::const_iterator d_it = den.begin(); d_it != den.end(); ++d_it)
      {

        // TODO check zero
        ss << "+("
           << d_it->coef.numerator.get_str() <<  "/"
           << d_it->coef.denominator.get_str() << ")";
        if (d_it->powers.front() > 0)
          ss << "*d^" << d_it->powers.front();

      }
  ss << ")";  
  return ss.str();
}


void Permutations::form_export(std::ostream & sout, const NumeratorTensor& ntens, bool compact)
{
  sout << "* " << endl;
  sout << "* " << endl;
  sout << "* " << ntens.info() << endl;
  sout << "* " << endl;
  sout << "* " << endl;

  sout << "*--#[ declare :" << endl;
  // loop configuration to b symbol
  map<LoopInfo, string> contr2b;
  // b symbol abreviations
  vector<string> b_abrev;

  auto rp_it =min_perms.begin();
  for (size_t i = 0; i < min_perms.size(); i++)
    {
      ostringstream bstr;
      bstr << "ttrB" << i+1;
      sout << "S " << bstr.str() << ";" << endl;
      b_abrev.push_back(bstr.str());
      contr2b[rp_it->first] = bstr.str();
      ++rp_it;
    }
  sout << "*--#] declare :" << endl;

  sout << "*--#[ reduce :" << endl;
  sout << " id " << ntens.id_lhs() << " = " << endl << endl;

  if (compact)
    {
      // Compact otput, may be memory consuming
      for (auto p_sigma: perms)
        {
          sout << " + " << ntens.g_perm(p_sigma) << " * (" << endl;

          // (sp)*(b) -> multiplicity
          unordered_map<string, size_t> sp_b_mult;
          for (auto p_tau: perms)
            {
              string sp_b = ntens.contract(p_tau) + " * " + contr2b[contract(p_sigma, p_tau)];
              
              if (auto k = sp_b_mult.find(sp_b); k != sp_b_mult.end())
                k->second += 1;
              else
                sp_b_mult[sp_b] = 1;
            }
          
          for( const auto& [sp_x_b, mul_factor] : sp_b_mult )
            sout << "\t+ " << sp_x_b << " * " << mul_factor << endl;
          
          sout << "   )" << endl;
          sout << endl;
        }
      
    }
  else
    {
      // Simple output
      // Loop over tensor structures
      for (auto p_sigma: perms)
        {
          // Loop over projector elements
          for (auto p_tau: perms)
            sout << " + " << ntens.g_perm(p_sigma)
                 << " * " << ntens.contract(p_tau)
                 << " * " << contr2b[contract(p_sigma, p_tau)]
                 << endl;
          sout << endl;
        }
    }
  sout << ";" << endl;
  sout << "*--#] reduce :" << endl;
 


  sout << "*--#[ subs :" << endl;
  // TODO
  
  for (size_t i = 0; i < b.size(); i++)
    {
      sout << " id " << b_abrev[i] << " = " << nd2rat(b[i]);
      if (min_multiplicity[i] > 1) sout << "/" << min_multiplicity[i];
      sout << ";" << endl;
    }
  sout << "*--#] subs :" << endl;
  
}



void usage()
{

  cout << "USAGE" << endl;
  cout << "\techo \"p1(m1)*p2(mu2)...*p2N(mu2N)\" | ./TTR [-qc] [-w workers] [-o outfile]" << endl;
  cout << endl;
  cout << "NAME" << endl;
  cout << "\t TTR = [T]adpoles [T]ensor [R]eduction" << endl;
  cout << endl;
  cout << "OPTIONS" << endl;
  cout << "  -q                         do not produce noisy otput during matrix inversion." << endl;
  cout << "  -c                         combine common output expressions to make it more compact." << endl;
  cout << "  -w WORKERS                 number of parallel jobs for matrix inversion, default=1." << endl;
  cout << "  -o OUTFILE                 output file name, if not specified used standard output." << endl;
  cout << endl;
  cout << "AUTHORS" << endl;
  cout << "\tAndrey Pikelner" << endl;
}



int main(int argc, char *argv[])
{
  try {

    size_t workers    = 1;
    bool out_to_file  = false;
    string ofname;
    bool need_compact = false;
    bool be_quiet     = false;

    for(;;)
      {
        switch(getopt(argc, argv, "o:w:cqh")) 
          {
          case 'o':
            out_to_file = true;
            ofname = optarg;
            continue;

          case 'w':
            workers = atoi(optarg);
            continue;

          case 'c':
            need_compact = true;
            continue;

          case 'q':
            be_quiet = true;
            continue;


          case 'h':
          default :
            usage();
            return 1;
            break;

          case -1:
            break;
          }

        break;
      }
    
    // read standard input
    string input;
    getline(cin, input);

    // Manipulations with input tensor
    NumeratorTensor ntens(input);

    // Representative permutations construction
    Permutations perms(ntens.rank(), workers, be_quiet);

    if (out_to_file)
      {
        cout << " - Saving result to \"" << ofname << "\"" << endl;
        std::ofstream fout(ofname);
        perms.form_export(fout, ntens, need_compact);
      }
    else
      {
        cout << " - Output result to standard output" << endl;
        perms.form_export(cout, ntens, need_compact);
      }

    cout << "All done, exiting!" << endl;
  }
  catch (std::exception &e)
    {
      std::cout << " ~ Exception : " << e.what() << std::endl;
    }
  
  return 0;
}
