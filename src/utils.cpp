#include "utils.h"

char*
re2post(char *re)
{
	int nalt, natom;
	static char buf[8000];
	char *dst;
	struct {
		int nalt;
		int natom;
	} paren[100], *p;
	
	p = paren;
	dst = buf;
	nalt = 0;
	natom = 0;
	if(strlen(re) >= sizeof buf/2)
		return NULL;
	for(; *re; re++){
		switch(*re){
		case '(':
			if(natom > 1){
				--natom;
				*dst++ = '.';
			}
			if(p >= paren+100)
				return NULL;
			p->nalt = nalt;
			p->natom = natom;
			p++;
			nalt = 0;
			natom = 0;
			break;
		case '|':
			if(natom == 0)
				return NULL;
			while(--natom > 0)
				*dst++ = '.';
			nalt++;
			break;
		case ')':
			if(p == paren)
				return NULL;
			if(natom == 0)
				return NULL;
			while(--natom > 0)
				*dst++ = '.';
			for(; nalt > 0; nalt--)
				*dst++ = '|';
			--p;
			nalt = p->nalt;
			natom = p->natom;
			natom++;
			break;
		case '*':
		case '+':
		case '?':
			if(natom == 0)
				return NULL;
			*dst++ = *re;
			break;
		default:
			if(natom > 1){
				--natom;
				*dst++ = '.';
			}
			*dst++ = *re;
			natom++;
			break;
		}
	}
	if(p != paren)
		return NULL;
	while(--natom > 0)
		*dst++ = '.';
	for(; nalt > 0; nalt--)
		*dst++ = '|';
	*dst = 0;
	return buf;
}


std::string stream_as_string(const std::string& path)
{
    std::filebuf fb;
    fb.open (path,std::ios::in);
    std::istream stm(&fb);
    std::string str;
    char c ;
    while( stm.get(c) ) str += c ;
    str.erase(std::remove(str.begin(), str.end(), '\n'), str.cend()); 
    return str;
}

int matches(const std::string& bin, std::regex reg, std::fstream& writefile)
{
    std::string match;
    std::sregex_iterator currentMatch(bin.begin(), bin.end(), reg);
    std::sregex_iterator last_match;
	int out = 0;
    while(currentMatch != last_match)
    {
        std::smatch match = *currentMatch;
        writefile<<match.str()<<"   "<<match.position()<<"\n";
        currentMatch++;
		out++;
    }
    return out;
}

std::string translate(const std::string& str)
{
  char * cstr = new char [str.length()+1];
  std::strcpy (cstr, str.c_str());
  std::string out = re2post(cstr);
  delete [] cstr;
  return out;
}

std::vector<char> getAlphabet(const std::string& regex)
{
  std::vector<char> alphabet{};
  for(auto e : regex)
  {
    if(e != '.' && e != '|' && e != '?' && e != '+' && e != '*')
    {
      alphabet.push_back(e);
    }
  }
  std::sort(alphabet.begin(), alphabet.end());
  auto last = std::unique(alphabet.begin(), alphabet.end());
  alphabet.erase(last, alphabet.end());
  return alphabet;
}


void matrixTotxt(const std::vector<std::vector<std::string>>& matrix, std::string& filename)
{
  std::fstream f;
  
  f.open(filename += ".txt", std::ios::out);
  if(f.good())
  {
    for(auto z : matrix)
    {
      for(auto s : z)
      {
        f << s <<" ";
      }
      f << "\n";
    }
  }
}
