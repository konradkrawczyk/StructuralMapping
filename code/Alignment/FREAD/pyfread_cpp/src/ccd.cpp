// Cyclic Coordinate Descent
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
using namespace std;


const float PI = 3.14159265358979323846f;


typedef class
{
  public:
    double xyz[3];
    double& operator[](int i);
} Atom;

double& Atom::operator[](int i)
{
  return xyz[i];
}

typedef double RotMat[3][3];


double unitX[3] = {1, 0, 0};
double unitZ[3] = {0, 0, 1};



//####################################################


vector<Atom> readCoord(FILE* file)
{
    if (!file)
    {
      fprintf(stderr, "ERROR reading input stream\n");
      exit(3);
    }
    
    vector<Atom> coord;
    char line[100];
    double x, y, z;
    
    while (fgets(line, 100, file))
    {
      if (!sscanf(line, "%lf %lf %lf", &x, &y, &z))
        break;
      Atom atm = {{x,y,z}};
      coord.push_back(atm);
    }
    
    return coord;
}

vector<Atom> readCoordFile(char * filename)
{
    vector<Atom> coord;
    char line[100];
    double x, y, z;

    FILE* file = fopen(filename, "r");

    if (!file)
    {
      fprintf(stderr, "ERROR reading file: %s\n", filename);
      exit(2);
    }
    
    while (fgets(line, 100, file))
    {
      if (!sscanf(line, "%lf %lf %lf", &x, &y, &z))
        break;
      Atom atm = {{x,y,z}};
      coord.push_back(atm);
    }

    fclose(file);
    
    return coord;
}

void printCoord(vector<Atom>& coord)
{
  for (unsigned i=0; i<coord.size(); ++i)
  {
    //printf("%10.3f %10.3f %10.3f\n", coord[i][0], coord[i][1], coord[i][2]);
    printf("%.4f\t%.4f\t%.4f\n", coord[i][0], coord[i][1], coord[i][2]);
  }
}


//####################################################


double sqrdist(double* a, double* b)
{
  return pow((a[0] - b[0]), 2) + pow((a[1] - b[1]), 2) + pow((a[2] - b[2]), 2);
}

double norm(double* a)
{
  return sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));
}

void cross(double* a, double* b, double* output)
{
  output[0] = a[1]*b[2] - a[2]*b[1];
  output[1] = a[2]*b[0] - a[0]*b[2];
  output[2] = a[0]*b[1] - a[1]*b[0];
}

double inner_product(double* a, double* b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


void assign(double* to, double* from)
{
  for (unsigned i=0; i<3; ++i)
    to[i] = from[i];
}

void subtract(double* a, double* b)
{
  a[0] -= b[0];
  a[1] -= b[1];
  a[2] -= b[2];
}

void divide(double* a, double x)
{
  a[0] /= x;
  a[1] /= x;
  a[2] /= x;
}

void rotate(double* a, RotMat& r)
{
  double temp[3];
  for (unsigned i=0; i<3; ++i)
  {
    temp[i] = inner_product(r[i], a);
  }
  a[0] = temp[0];
  a[1] = temp[1];
  a[2] = temp[2];
}


//####################################################


double TerminalRMSD(vector<Atom>& loop, vector<Atom>& terminal)
{
    int loop_offset = loop.size()-terminal.size();
    double total = 0.0;
    for (unsigned i=0; i<terminal.size(); ++i)
    {
      total += sqrdist(loop[loop_offset+i].xyz, terminal[i].xyz);
    }
    return sqrt(total / terminal.size());
}

void Translate(vector<Atom>& coords, double* T, int start)
{
    for (unsigned i=start; i<coords.size(); ++i)
      subtract(coords[i].xyz, T);
}

void SetRotMatY(double angle, RotMat& rotmat)
{
    rotmat[0][0] = cos(angle);
    rotmat[0][1] = 0;
    rotmat[0][2] = -sin(angle);
    rotmat[1][0] = 0;
    rotmat[1][1] = 1;
    rotmat[1][2] = 0;
    rotmat[2][0] = sin(angle);
    rotmat[2][1] = 0;
    rotmat[2][2] = cos(angle);
}
void SetRotMatZ(double angle, RotMat& rotmat)
{
    rotmat[0][0] = cos(angle);
    rotmat[0][1] = -sin(angle);
    rotmat[0][2] = 0;
    rotmat[1][0] = sin(angle);
    rotmat[1][1] = cos(angle);
    rotmat[1][2] = 0;
    rotmat[2][0] = 0;
    rotmat[2][1] = 0;
    rotmat[2][2] = 1;
}

void Rot(RotMat& rotmat, vector<Atom>& coords, int pos)
{
    rotate(coords[pos].xyz, rotmat);
}

void ListRot(RotMat& rotmat, vector<Atom>& coords, int start)
{
    for (unsigned i=start; i<coords.size(); ++i)
    {
      rotate(coords[i].xyz, rotmat);
    }
}

double OptimalAngle(vector<Atom>& move, vector<Atom>& fix)
{
    int move_offset = move.size()-3;
    double numerator = 0;
    double denominator = 0;
    for (int i=0; i<3; ++i)
    {
      double* f_vec = fix[i].xyz;
      
      double r_norm = norm(move[move_offset+i].xyz);
      
      double r[3] = {move[move_offset+i][0], move[move_offset+i][1], 0};
      divide(r, norm(r));
      
      double s[3];
      cross(r, unitZ, s);
      
      numerator += inner_product(s, f_vec)*r_norm;
      denominator += inner_product(r, f_vec)*r_norm;
    }
    return atan(numerator/denominator);
}



//####################################################



void LoopClosure(vector<Atom>& loop, int n, double* origin, vector<Atom>& terminal)
{   
    Translate(loop, origin, n);
    Translate(terminal, origin, 0);
    
    
    //# #############
    //# Align the second atom on the X-Z plane
    //# #############
    
    double axis[3];
    assign(axis, loop[n+1].xyz);
    subtract(axis, loop[n].xyz);
    axis[2] -= loop[n+1].xyz[2];
    divide(axis, norm(axis));
    
    double phi = atan2(axis[1], axis[0]);
    
    RotMat rotmatZ;
    SetRotMatZ(-phi, rotmatZ);
    
    Rot(rotmatZ, loop, n+1);
    ListRot(rotmatZ, loop, n+2);
    ListRot(rotmatZ, terminal, 0);
    
    //# #############
  
  
    //# #############
    //# Align the second atom on Z-axis
    //# #############
    
    assign(axis, loop[n+1].xyz);
    subtract(axis, loop[n].xyz);
    divide(axis, norm(axis));
    
    double theta = (PI / 2) + atan2(axis[2], axis[0]);
    
    RotMat rotmatY;
    SetRotMatY(-theta, rotmatY);

    Rot(rotmatY, loop, n+1);
    ListRot(rotmatY, loop, n+2);
    ListRot(rotmatY, terminal, 0);
    
    //# #############
    
    
    double psi = OptimalAngle(loop, terminal);
    RotMat rotmatZ2;
    SetRotMatZ(-psi, rotmatZ2);
    ListRot(rotmatZ2, loop, n+1);
    
    
    //# ###################
    //# Putting coordinates to the original positions
    //# ###################
    SetRotMatZ(phi, rotmatZ);
    SetRotMatY(theta, rotmatY);
    
    Rot(rotmatY, loop, n+1);
    ListRot(rotmatY, loop, n+2);
    ListRot(rotmatY, terminal, 0);
    Rot(rotmatZ, loop, n+1);
    ListRot(rotmatZ, loop, n+2);
    ListRot(rotmatZ, terminal, 0);
    divide(origin, -1); // CAUTION, THIS HAS SIDE EFFECTS
    Translate(loop, origin, n);
    Translate(terminal, origin, 0);
}


void CyclicCoordinateDescent(vector<Atom>& loop, vector<Atom>& terminal, int max_iterations, double max_rmsd, double& out_rmsd, int& out_iterations)
{
    double rmsd = 99999999;
    int iter=0;
    bool closed = false;
    while (!closed && (iter < max_iterations))
    {
      for (unsigned i=3; i<loop.size()-3; ++i)
      {
        if ((i-3)%3 == 2)
        {
          // omega angle
        }
        else
        {
          double origin[3];
          assign(origin, loop[i].xyz);
          LoopClosure(loop, i, origin, terminal);
          rmsd = TerminalRMSD(loop, terminal);
          if (rmsd <= max_rmsd)
          {
            closed = true;
            break;
          }
        }
      }
      iter += 1;
      if (rmsd != rmsd)
      {
        fprintf(stderr, "RMSD == NaN!\n");
        break;
      }
    }
    out_rmsd = rmsd;
    out_iterations += iter;
}


//####################################################


void reverse(vector<Atom>& atoms)
{
  double temp[3];
  int size=atoms.size();
  for (int i=0; i<size/2; ++i)
  {
    int j = size - 1 - i;
    assign(temp, atoms[j].xyz);
    assign(atoms[j].xyz, atoms[i].xyz);
    assign(atoms[i].xyz, temp);
  }
}


void close_loop(vector<Atom>& Nterm, vector<Atom>& loop, vector<Atom>& Cterm, int max_iterations, double max_rmsd, double& rmsd, int& iterations)
{
    const int start_iter = 1;
    double n_rmsd = -1;
    double c_rmsd = -1;
    int n_iterations = start_iter;
    int c_iterations = start_iter;
    
    bool nclosed = false;
    bool cclosed = false;
    bool loop_reversed = false;
    
    reverse(Nterm);
    
    int remaining_iterations = max_iterations + 2*start_iter;
    while (!(nclosed && cclosed) && (remaining_iterations > 0))
    {
//       int maxstepiter = remaining_iterations;
//       if (!nclosed && !cclosed)
//         maxstepiter = (maxstepiter + 1) / 2;
      
      if (!cclosed)
      {
        if (loop_reversed)
        {
          reverse(loop);
          loop_reversed = !loop_reversed;
        }
        // Close C-terminal
        // CyclicCoordinateDescent(loop, Cterm, min(c_iterations, maxstepiter), max_rmsd, c_rmsd, c_iterations);
        CyclicCoordinateDescent(loop, Cterm, c_iterations, max_rmsd, c_rmsd, c_iterations);
        cclosed = (c_rmsd <= max_rmsd);
      }
      
      if (!nclosed)
      {
        if (!loop_reversed)
        {
          reverse(loop);
          loop_reversed = !loop_reversed;
        }
        // CyclicCoordinateDescent(loop, Nterm, min(n_iterations, maxstepiter), max_rmsd, n_rmsd, n_iterations);
        CyclicCoordinateDescent(loop, Nterm, n_iterations, max_rmsd, n_rmsd, n_iterations);
        nclosed = (n_rmsd <= max_rmsd);
      }
      
      remaining_iterations = max_iterations + 2*start_iter - n_iterations - c_iterations;
    }
    
    if (loop_reversed)
    {
      reverse(loop);
      loop_reversed = !loop_reversed;
    }
    
    reverse(Nterm);
    
    // Combined RMSD and iterations for N and C termini
    rmsd = sqrt((pow(n_rmsd, 2) + pow(c_rmsd, 2)) / 2.0);
    iterations = n_iterations + c_iterations - 2*start_iter;
}


//####################################################


int main(int argc, char** argv)
{
  if ((argc != 6) && (argc != 3))
  {
      fprintf(stderr, "USAGE: ccd <target_anchor_rmsd> <max_iterations> [<N-anchor> <loop> <C-anchor>]\n");
      exit(2);
  }
  
  double max_rmsd = atof(argv[1]);
  int max_iterations = atoi(argv[2]);
  
  vector<Atom> Nanchor;
  vector<Atom> loop;
  vector<Atom> Canchor;
  
  
  if (argc == 3)
  {
    Nanchor = readCoord(stdin);
    if (Nanchor.empty())
    {
      fprintf(stderr, "ERROR: No N anchor!\nUSAGE: ccd [<N-anchor> <loop> <C-anchor>]\n");
      exit(3);
    }
    loop = readCoord(stdin);
    if (loop.empty())
    {
      fprintf(stderr, "ERROR: No loop coordinates!\nUSAGE: ccd [<N-anchor> <loop> <C-anchor>]\n");
      exit(3);
    }
    Canchor = readCoord(stdin);
    if (Canchor.empty())
    {
      fprintf(stderr, "ERROR: No C anchor!\nUSAGE: ccd [<N-anchor> <loop> <C-anchor>]\n");
      exit(3);
    }
  }
  else
  {
    Nanchor = readCoordFile(argv[3]);
    loop = readCoordFile(argv[4]);
    Canchor = readCoordFile(argv[5]);
  }
  
  
  if (Nanchor.size() != Canchor.size())
  {
    fprintf(stderr, "ERROR: Anchor lengths not equal!\n");
    exit(4);
  }
  if (Nanchor.size() != 3)
  {
    fprintf(stderr, "ERROR: Anchors need to have 3 atoms (N, CA, C)!\n");
    exit(5);
  }
  if (loop.size() == 0)
  {
    fprintf(stderr, "ERROR: No loop coordinates!\n");
    exit(6);
  }
  if (loop.size() % 3 != 0)
  {
    fprintf(stderr, "ERROR: Loop needs to have 3 atoms per residue (N, CA, C)!\n");
    exit(6);
  }
  
  
  double rmsd = -1;
  int iterations = 0;
  
  close_loop(Nanchor, loop, Canchor, max_iterations, max_rmsd, rmsd, iterations);
  
  printCoord(loop);
  
  fprintf(stderr, "%d\n", iterations);
  fprintf(stderr, "%.4f\n", rmsd);
  
  
  return 1 - (rmsd >= 0 && rmsd <= max_rmsd);
}
