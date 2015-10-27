#include "set_parameters.h"

Parameters::Parameters()
{

    num_nodes=unlikely;
    average_k=unlikely;
    max_degree=unlikely;

    tau=2;
    tau2=1;

    mixing_parameter=unlikely;
    mixing_parameter2=unlikely;

    beta=1.5;

    overlapping_nodes=0;
    overlap_membership=0;

    nmin=unlikely;
    nmax=unlikely;

    randomf=false;
    fixed_range=false;
    excess=false;
    defect=false;

    clustering_coeff=unlikely;

    command_flags.push_back("-N");			//0
    command_flags.push_back("-k");			//1
    command_flags.push_back("-maxk");		//2
    command_flags.push_back("-mut");		//3
    command_flags.push_back("-t1");			//4
    command_flags.push_back("-t2");			//5
    command_flags.push_back("-minc");		//6
    command_flags.push_back("-maxc");		//7
    command_flags.push_back("-on");			//8
    command_flags.push_back("-om");			//9
    command_flags.push_back("-beta");		//10
    command_flags.push_back("-muw");		//11
    command_flags.push_back("-C");			//12

}

void Parameters::set_random()
{

    FILE_LOG(logINFO)<<"this is a random network";
    mixing_parameter=0;
    mixing_parameter2=0;
    overlapping_nodes=0;
    overlap_membership=0;
    nmax=num_nodes;
    nmin=num_nodes;

    fixed_range=true;
    excess=false;
    defect=false;
}

bool Parameters::arrange()
{

    if(randomf)
        set_random();

    if (num_nodes==unlikely)
    {
        cerr<<"\n***********************\nERROR:\t number of nodes unspecified";
        return false;
    }

    if (average_k==unlikely)
    {
        cerr<<"\n***********************\nERROR:\t average degree unspecified";
        return false;
    }

    if (max_degree==unlikely)
    {
        cerr<<"\n***********************\nERROR:\t maximum degree unspecified";
        return false;
    }

    if (mixing_parameter2==unlikely)
    {
        cerr<<"\n***********************\nERROR:\t weight mixing parameter (option -muw) unspecified";
        return false;
    }

    if(mixing_parameter==unlikely)
        mixing_parameter=mixing_parameter2;

    if(overlapping_nodes<0 || overlap_membership<0)
    {
        FILE_LOG(logERROR)<<"ERROR:\tsome positive parameters are negative";
        return -1;
    }

    if (num_nodes<=0 || average_k<=0 || max_degree<=0 || mixing_parameter<0 || mixing_parameter2<0 || (nmax<=0 && nmax!=unlikely) || (nmin<=0 && nmin!=unlikely) )
    {
        FILE_LOG(logERROR)<<"ERROR:\tsome positive parameters are negative";
        return -1;
    }
    if(mixing_parameter > 1 || mixing_parameter2 > 1)
    {
        FILE_LOG(logERROR)<<"ERROR:\tmixing parameter > 1 (must be between 0 and 1)";
        return -1;
    }

    if(nmax!= unlikely && nmin!=unlikely)
        fixed_range=true;
    else
        fixed_range=false;

    if(excess && defect)
    {
        FILE_LOG(logERROR)<<"ERROR:\tboth options -inf and -sup cannot be used at the same time";
        return false;
    }

    FILE_LOG(logINFO)<<"number of nodes:\t"<<num_nodes;
    FILE_LOG(logINFO)<<"average degree:\t"<<average_k;
    FILE_LOG(logINFO)<<"maximum degree:\t"<<max_degree;
    FILE_LOG(logINFO)<<"exponent for the degree distribution:\t"<<tau;
    FILE_LOG(logINFO)<<"exponent for the community size distribution:\t"<<tau2;
    FILE_LOG(logINFO)<<"mixing parameter(topology):\t"<<mixing_parameter;
    FILE_LOG(logINFO)<<"mixing parameter (weights):\t"<<mixing_parameter2;
    FILE_LOG(logINFO)<<"beta exponent:\t"<<beta;
    FILE_LOG(logINFO)<<"number of overlapping nodes:\t"<<overlapping_nodes;
    FILE_LOG(logINFO)<<"number of memberships of the overlapping nodes:\t"<<overlap_membership;
    if(clustering_coeff!=unlikely)
        FILE_LOG(logINFO) << "Average clustering coefficient: "<<clustering_coeff;

    if (fixed_range)
    {
        FILE_LOG(logINFO)<<"Community size range set equal to ["<<nmin<<" , "<<nmax<<"]";

        if (nmin>nmax)
        {
            FILE_LOG(logERROR)<<"ERROR: INVERTED COMMUNITY SIZE BOUNDS";
            return false;
        }

        if(nmax>num_nodes)
        {
            FILE_LOG(logERROR)<<"ERROR: maxc BIGGER THAN THE NUMBER OF NODES";
            return false;
        }
    }

    return true;

}

/**
 * @brief Parameters::set
 * @param flag
 * @param num
 * @return
 */
bool Parameters::set(string & flag, string & num)
{
    // false is something goes wrong
    FILE_LOG(logINFO)<<"setting... "<<flag<<" "<<num;
    double err;
    if(!cast_string_to_double(num, err))
    {
        cerr<<"\n***********************\nERROR while reading parameters";
        return false;
    }

    if (flag==command_flags[0])
    {
        if (fabs(err-int (err))>1e-8)
        {
            cerr<<"\n***********************\nERROR: number of nodes must be an integer";
            return false;
        }
        num_nodes=cast_int(err);
    }
    else if(flag==command_flags[1])
    {
        average_k=err;
    }
    else if(flag==command_flags[2])
    {
        max_degree=cast_int(err);
    }
    else if(flag==command_flags[3])
    {
        mixing_parameter=err;
    }
    else if(flag==command_flags[11])
    {
        mixing_parameter2=err;
    }
    else if(flag==command_flags[10])
    {
        beta=err;
    }
    else if(flag==command_flags[4])
    {
        tau=err;
    }
    else if(flag==command_flags[5])
    {
        tau2=err;
    }

    else if(flag==command_flags[6])
    {
        if (fabs(err-int (err))>1e-8)
        {
            cerr<<"\n***********************\nERROR: the minumum community size must be an integer";
            return false;
        }
        nmin=cast_int(err);
    }
    else if(flag==command_flags[7])
    {
        if (fabs(err-int (err))>1e-8)
        {
            cerr<<"\n***********************\nERROR: the maximum community size must be an integer";
            return false;
        }
        nmax=cast_int(err);
    }
    else if(flag==command_flags[8])
    {
        if (fabs(err-int (err))>1e-8)
        {
            cerr<<"\n***********************\nERROR: the number of overlapping nodes must be an integer";
            return false;
        }
        overlapping_nodes=cast_int(err);
    }
    else if(flag==command_flags[9])
    {
        if (fabs(err-int (err))>1e-8)
        {

            cerr<<"\n***********************\nERROR: the number of membership of the overlapping nodes must be an integer";
            return false;
        }
        overlap_membership=cast_int(err);
    }
    else if(flag==command_flags[12])
    {
        clustering_coeff=err;
    }
    else
    {
        cerr<<"\n***********************\nERROR while reading parameters: "<<flag<<" is an unknown option";
        return false;
    }
    return true;
}

/**
 * @brief print_usage
 */
void print_usage()
{

    FILE_LOG(logINFO)<<"\nTo run the program type \n./benchmark [FLAG] [P]";
    FILE_LOG(logINFO)<<"\n----------------------\n";
    FILE_LOG(logINFO)<<"To set the parameters, type:"<<endl;
    FILE_LOG(logINFO)<<"-N\t\t[number of nodes]";
    FILE_LOG(logINFO)<<"-k\t\t[average degree]";
    FILE_LOG(logINFO)<<"-maxk\t\t[maximum degree]";
    FILE_LOG(logINFO)<<"-mut\t\t[mixing parameter for the topology]";
    FILE_LOG(logINFO)<<"-muw\t\t[mixing parameter for the weights]";
    FILE_LOG(logINFO)<<"-beta\t\t[exponent for the weight distribution]";
    FILE_LOG(logINFO)<<"-t1\t\t[minus exponent for the degree sequence]";
    FILE_LOG(logINFO)<<"-t2\t\t[minus exponent for the community size distribution]";
    FILE_LOG(logINFO)<<"-minc\t\t[minimum for the community sizes]";
    FILE_LOG(logINFO)<<"-maxc\t\t[maximum for the community sizes]";
    FILE_LOG(logINFO)<<"-on\t\t[number of overlapping nodes]";
    FILE_LOG(logINFO)<<"-om\t\t[number of memberships of the overlapping nodes]";
    FILE_LOG(logINFO)<<"-C\t\t[Average clustering coefficient]";
    FILE_LOG(logINFO)<<"----------------------\n";
    FILE_LOG(logINFO)<<"It is also possible to set the parameters writing flags and relative numbers in a file. To specify the file, use the option:";
    FILE_LOG(logINFO)<<"-f\t[filename]";
    FILE_LOG(logINFO)<<"You can set the parameters both writing some of them in the file, and using flags from the command line for others."<<endl;
    FILE_LOG(logINFO)<<"-N, -k, -maxk, -muw have to be specified. For the others, the program can use default values:";
    FILE_LOG(logINFO)<<"t1=2, t2=1, on=0, om=0, beta=1.5, mut=muw, minc and maxc will be chosen close to the degree sequence extremes.";
    FILE_LOG(logINFO)<<"If you don't specify -C the rewiring process for raising the average clustering coefficient will not be performed";
    FILE_LOG(logINFO)<<"If you set a parameter twice, the latter one will be taken.";
    FILE_LOG(logINFO)<<"\n-------------------- Other options ---------------------------\n";
    FILE_LOG(logINFO)<<"To have a random network use:";
    FILE_LOG(logINFO)<<"-rand";
    FILE_LOG(logINFO)<<"Using this option will set muw=0, mut=0, and minc=maxc=N, i.e. there will be one only community.";
    FILE_LOG(logINFO)<<"Use option -sup (-inf) if you want to produce a benchmark whose distribution of the ratio of external degree/total degree ";
    FILE_LOG(logINFO)<<"is superiorly (inferiorly) bounded by the mixing parameter.";
    FILE_LOG(logINFO)<<"\n-------------------- Examples ---------------------------\n";
    FILE_LOG(logINFO)<<"Example1:";
    FILE_LOG(logINFO)<<"./benchmark -N 1000 -k 15 -maxk 50 -muw 0.1 -minc 20 -maxc 50";
    FILE_LOG(logINFO)<<"Example2:";
    FILE_LOG(logINFO)<<"./benchmark -f flags.dat -t1 3";
    FILE_LOG(logINFO)<<"\n-------------------- Other info ---------------------------\n";
    FILE_LOG(logINFO)<<"Read file ReadMe.txt for more info."<<endl;
}

bool set_from_file(string & file_name, Parameters & par1)
{
    int h= file_name.size();
    char b[h+1];
    cast_string_to_char(file_name, b);

    ifstream in(b);
    if (!in.is_open())
    {
        FILE_LOG(logERROR)<<"File "<<file_name<<" not found. Where is it?";
        return false;
    }

    string temp;
    while(in>>temp)  			// input file name
    {
        if(temp=="-rand")
            par1.randomf=true;
        else if(temp=="-sup")
            par1.excess=true;
        else if(temp=="-inf")
            par1.defect=true;
        else
        {
            string temp2;
            in>>temp2;
            if(temp2.size()>0)
            {
                if(temp=="-f" && temp2!=file_name)
                {
                    if(set_from_file(temp2, par1)==false)
                        return false;
                }
                if(temp!="-f")
                {
                    if(par1.set(temp, temp2)==false)
                        return false;
                }
            }
            else
            {
                FILE_LOG(logERROR)<<"Errow while reading parameters";
                return false;
            }
        }
    }
    return true;
}

bool set_parameters(int argc, char * argv[], Parameters & par1)
{
    if (argc <= 1)   // if no arguments, return statement about program usage.
    {
        print_usage();
        return false;
    }

    int argct = 0;
    string temp;
    while (++argct < argc)  			// input file name
    {
        temp = argv[argct];
        if(temp=="-rand")
            par1.randomf=true;
        else if(temp=="-sup")
            par1.excess=true;
        else if(temp=="-inf")
            par1.defect=true;
        else
        {
            argct++;
            string temp2;
            if(argct<argc)
            {
                temp2 = argv[argct];
                if(temp=="-f")
                {
                    if(set_from_file(temp2, par1)==false)
                        return false;
                }
                if(temp!="-f")
                {
                    if(par1.set(temp, temp2)==false)
                        return false;
                }
            }
            else
            {
                FILE_LOG(logERROR)<<"ERROR while reading parameters";
                return false;
            }
        }
    }

    if(par1.arrange()==false)
        return false;

    return true;
}

