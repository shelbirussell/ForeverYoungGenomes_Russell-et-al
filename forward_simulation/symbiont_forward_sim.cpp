#include <string>
#include <string.h>
#include <iostream>
#include <time.h> 
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <map>
#include <set>
#include <utility>
#include <list>
#include <ctype.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

// To check for gsl dependencies in different architectures: gsl-config --libs
// So, on Odyssey, the compile command is: g++ -std=c++0x -O3 -lgsl -lgslcblas -lm -o SYM_REC SYM_REC.cpp

// g++ -std=c++0x -O3 -lgsl -ltcmalloc -o SYM_REC SYM_REC.cpp

using namespace std;

/// this is the declaration of the random number generator to be used throughout
const gsl_rng *rng;

/// command line information and global parameters
class cmd_line {
public:
    
    /// generations to run for
    int generations ;
    
    /// diploid effective population size
    int host_ne ;
    
    /// growth rate
    long double growth_rate ;
    
    /// haploid effective population sizes through time
    vector<int> symbiont_ne ;

    /// mutation rate per symbiont genome
    long double mu ;
    
    /// mean length of exponentially distributed recombination tracts
    long double r ;
    long double rho ;
    
    /// seed for rng
    long double seed ;
    
    /// horizontal transmission rate per symbiont lineage
    long double horizontal ;
    
    /// print frequency
    int print_freq ;
    
    /// read relevant information
    void read_cmd_line ( int argc, char *argv[] ) ;
    
} ;

void compute_symbiont_ne ( long double n0, long double rate, vector<int> &symbiont_ne ) {
    for ( int i = 0 ; i <= 3 ; i ++ ) {
        symbiont_ne.push_back( n0*pow(rate,i) ) ;
    }
}

void cmd_line::read_cmd_line ( int argc, char *argv[] ) {
    
    host_ne = 500 ;
    growth_rate = 2 ;
    compute_symbiont_ne( 100, growth_rate, symbiont_ne ) ;
    mu = 0.02 ;
    rho = 0.0 ;
    r = 0.0 ;
    generations = 10001 ;
    horizontal = 0.01 ;
    print_freq = 1000 ;
    
    /// accept command line parameters
    for (int i=1; i<argc; i++) {
        if ( strcmp(argv[i],"-n") == 0 ) {
            host_ne = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-p") == 0 ) {
            print_freq = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-r") == 0 ) {
            growth_rate = atof(argv[++i]) ;
            symbiont_ne.clear() ;
            compute_symbiont_ne(100,growth_rate,symbiont_ne) ;
        }
        if ( strcmp(argv[i],"-g") == 0 ) {
            generations = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-u") == 0 ) {
            mu = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-l") == 0 ) {
            r = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--rho") == 0 ) {
            rho = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-s") == 0 ) {
            seed = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-h") == 0 ) {
            horizontal = atof(argv[++i]) ;
        }
    }
}

void create_population( vector<vector<vector<float> *> > &host, vector<float> *individual, int ne ) {
    vector<vector<float> *> host_individual (100,individual) ;
    for ( int i = 0 ; i < ne ; i ++ ) {
        host.push_back(host_individual) ;
    }
}

void clone_and_mutate ( vector<vector<float> *> &host, int ne, cmd_line &options, list<vector<float> *> &all_inds ) {
    
    /// first clone all individuals
    vector<vector<float> *> new_host ;
    for ( int i = 0 ; i < ne ; i ++ ) {
        new_host.push_back(host.at(gsl_rng_uniform_int(rng, host.size()))) ;
    }
    
    /// now mutate as needed
    for ( int i = 0 ; i < new_host.size() ; i ++ ) {
        int mutations = gsl_ran_poisson( rng, options.mu ) ;
        if ( mutations > 0 ) {
            vector<float> *symbiont = new vector<float> () ;
            *symbiont = *new_host[i] ;
            for ( int m = 0 ; m < mutations ; m ++ ) {
                symbiont->push_back(gsl_rng_uniform(rng)) ;
            }
            sort(symbiont->begin(),symbiont->end()) ;
            all_inds.push_back( symbiont ) ;
            swap( new_host[i], symbiont ) ;
        }
    }
    
    /// swap pointers
    swap(host,new_host) ;
}

void migrate ( vector<vector<vector<float> *> > &host, cmd_line &options ) {
    for ( int h = 0 ; h < host.size() ; h ++ ) {
        for ( int s = 0 ; s < host[h].size() ; s ++ ) {
            if ( gsl_ran_bernoulli(rng,options.horizontal) == 1 ) {
                host[h][s] = host[gsl_rng_uniform_int(rng,host.size())][gsl_rng_uniform_int(rng,host[h].size())] ;
            }
        }
    }
}

void recombine( vector<vector<vector<float> *> > &host, cmd_line &options, list<vector<float> *> &all_inds ) {
    for ( int h = 0 ; h < host.size() ; h ++ ) {
        for ( int s = 0 ; s < host[h].size() ; s ++ ) {
            if ( gsl_ran_bernoulli(rng,options.rho) == 1 ) {
                
                /// donor information
                int d = gsl_rng_uniform_int(rng,host[h].size()) ;
                if ( *host[h][d] == *host[h][s] ) {
                    continue ;
                }
                
                /// tract information
                float length = gsl_ran_exponential( rng, options.r ) ;
                float start = gsl_rng_uniform( rng ) ;
                
                /// just replace lineage if the tract is greater than 1
                if ( length >= 1 ) {
                    host[h][s] = host[h][d] ;
                    continue ;
                }
                
                /// create new individual
                vector<float> *new_symbiont = new vector<float> () ;
                
                /// if the length is less than one but overlaps the edge
                if ( length + start > 1 ) {
                    for ( int m = 0 ; m < host[h][d]->size() ; m ++ ) {
                        if ( (*host[h][d])[m] < length + start - 1 ) {
                            new_symbiont->push_back( (*host[h][d])[m] ) ;
                        }
                        else {
                            break ;
                        }
                    }
                }
                /// now from beginning of existing to start
                for ( int m = 0 ; m < host[h][s]->size() ; m ++ ) {
                    if ( (*host[h][s])[m] > start ) {
                        break ;
                    }
                    if ( (*host[h][s])[m] > max(float(0),length+start-1) ) {
                        new_symbiont->push_back( (*host[h][s])[m] ) ;
                    }
                }
                /// now from start to end
                for ( int m = 0 ; m < host[h][d]->size() ; m ++ ) {
                    if ( (*host[h][d])[m] > start + length ) {
                        break ;
                    }
                    if ( (*host[h][d])[m] > start ) {
                        new_symbiont->push_back( (*host[h][d])[m] ) ;
                    }
                }
                /// finally the last strand of existing
                for ( int m = 0 ; m < host[h][s]->size() ; m ++ ) {
                    if ( (*host[h][s])[m] > length+start ) {
                        new_symbiont->push_back( (*host[h][s])[m] ) ;
                    }
                }
                all_inds.push_back( new_symbiont ) ;
                host[h][s] = new_symbiont ;
            }
        }
    }
}

void print_stats( vector<vector<vector<float> *> > &host, int g ) {
    
    map<float,vector<int> > allele_counts ;
    for ( int h = 0 ; h < host.size() ; h ++ ) {
        for ( int s = 0 ; s < 200 ; s ++ ) {
            for ( int m = 0 ; m < (*host[h][s]).size() ; m ++ ) {
                if ( allele_counts.find( (*host[h][s])[m] ) == allele_counts.end() ) {
                    allele_counts[(*host[h][s])[m]].resize(host.size(),0) ;
                }
                allele_counts[(*host[h][s])[m]][h] ++ ;
            }
        }
    }
    
	// print output for each allele:
    for ( auto m = allele_counts.begin() ; m != allele_counts.end() ; m ++ ) {
        cout << g << "\t" << m->first ;
        for ( int s = 0 ; s < m->second.size() ; s ++ ) {
            cout << "\t" << m->second.at(s) ;
        }
        cout << endl ;
    }
}

// select clams to reproduce via wright-fisher sampling
void select_hosts ( vector<vector<vector<float> * > > &hosts ) {
    
    vector<vector<vector<float> * > > new_hosts ;
    for ( int h = 0 ; h < hosts.size() ; h ++ ) {
        new_hosts.push_back( hosts[gsl_rng_uniform_int(rng, hosts.size())]) ;
    }
    
    swap( hosts, new_hosts ) ;
}

void collect_garbage ( list<vector<float> *> &all_inds, vector<vector<vector<float> * > > &host ) {
    
    // map to hold the individuals still in use
    map<vector<float>*, bool> extant ;
    for ( int h = 0 ; h < host.size() ; h ++ ) {
        for ( int s = 0 ; s < host[h].size() ; s ++ ) {
            extant[host[h][s]] = true ;
        }
    }
    
    for ( auto i = all_inds.begin() ; i != all_inds.end() ; ) {
        if ( !extant[*i] ) {
            delete *i ;
            i = all_inds.erase( i ) ;
        }
        else {
            i ++ ;
        }
    }
}

int main ( int argc, char **argv ) {
    
    //// read command line options
    cmd_line options ;
    options.read_cmd_line( argc, argv ) ;
    
    //initialize rng for gsl lookup table
    rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng, (long) options.seed);

    //create base individual
    vector<float> *individual = new vector<float>() ;
    list<vector<float> *> all_inds ;
    all_inds.push_back( individual ) ;
    
    //create population
    vector<vector<vector<float>* > > host ;
    create_population( host, individual, options.host_ne ) ;
    
    /// evolve in forward time
    for ( int g = 1 ; g < options.generations ; g ++ ) {

    	/// first select clams to reproduce via wright-fisher
        select_hosts ( host ) ;
                    
        /// for all steps within the generation
        for ( int step = 0 ; step < options.symbiont_ne.size() ; step ++ ) {
        
            /// first clone and mutate
            for ( int i = 0 ; i < host.size() ; i ++ ) {
                clone_and_mutate( host.at(i), options.symbiont_ne.at(step), options, all_inds ) ;
            }
            
            /// now migrate
            migrate( host, options ) ;
            
            /// now recombine
            recombine( host, options, all_inds ) ;
            
            /// garbage collect when pops are small
            if ( step == 0 ) {
                collect_garbage( all_inds, host ) ;
            }
        }
        
        // if the number of generations (g) divided by the print_freq equals zero, then print stats
        if ( g % options.print_freq == 0 ) {
            print_stats( host, g ) ;
        }
    }
    
    return(0) ;
}
