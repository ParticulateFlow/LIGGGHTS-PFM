#ifdef _WIN32
 #define _USE_MATH_DEFINES
 #include <math.h>
#endif

#include "liggghts_forcechainreader.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTable.h"
#include "vtkVariant.h"
#include "vtkObject.h"
#include "vtkLine.h"

#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include "vtkObjectFactory.h"
#include "vtkByteSwap.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkStringArray.h"
#include "vtkCharArray.h"
#include "vtkFloatArray.h"

#include <vtkCellData.h>

#include <algorithm>
#include <vector>
#include <string>
#include <sstream>

//#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"


//vtkCxxRevisionMacro(liggghts_forcechainreader, "$Revision: 2.0 $");
vtkStandardNewMacro(liggghts_forcechainreader);
//vtkInformationKeyMacro(liggghts_forcechainreader, TS_KEY, Integer);

//tiny little helper
void searchAndReplace(std::string& value, std::string const& search,std::string const& replace)
{
    std::string::size_type next;

    for(next = value.find(search);        // find first match
        next != std::string::npos;        // next is npos if search string was not found
        next = value.find(search,next))   // search for the next match starting after
                                          // the last match that was found.
    {
        // Inside the loop. So we found a match.
        value.replace(next,search.length(),replace);   // Do the replacement.
        next += replace.length();                      // Move to just after the replace
        // This is the point were we start
        // the next search from.
    }
}

liggghts_forcechainreader::liggghts_forcechainreader()
{
    this->FileName = 0;
    this->File = 0;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

liggghts_forcechainreader::~liggghts_forcechainreader()
{
    if (this->File)
    {
        this->File->close();
        delete this->File;
        this->File = NULL;
    }

    this->SetFileName(0);
    this->FileName = NULL;
}

void liggghts_forcechainreader::OpenFile()
{
    if (!this->FileName)
    {
        vtkErrorMacro(<<"FileName must be specified.");
        return;
    }
    // If the file was open close it.
    if (this->File)
    {
        this->File->close();
        delete this->File;
        this->File = NULL;
    }

    // Open the new file.

#ifdef _WIN32
    this->File = new ifstream(this->FileName, ios::in | ios::binary);
#else
    this->File = new ifstream(this->FileName, ios::in);
#endif
    if (! this->File || this->File->fail())
    {
        vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
        return;
    }
}

// remove special characters at the end
void trim(char *str)
{
    size_t i = strlen(str)-1;
    while( (i>=0) && ((str[i] == '\r') || (str[i] == '\n') || (str[i] == ' ')))
    {
        str[i] = '\0';
        i--;
    }
}

int liggghts_forcechainreader::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    vtkCellArray *lines;
    vtkPoints *points;
    vtkLine *tmpLine;

    this->OpenFile();
    char line[512]; // increase to get more chars !

    this->File->getline(line,sizeof(line)); // 1st line
    this->File->getline(line,sizeof(line)); // 2nd line = Timestep
    int TS=atoi(line);

    this->File->getline(line,sizeof(line)); // # of items
    this->File->getline(line,sizeof(line)); // 4th line = #Atoms
    int COUNT = atoi(line);
    //if (COUNT<1) return 1;

    double SHEAR=0;

    this->File->getline(line,sizeof(line)); //5th line = ITEM: BOX Bounds OR ITEM: ENTRIES
    trim(line);
    if (strncmp(line,"ITEM: BOX BOUNDS",15) == 0) {
        //override ..
        this->File->getline(line,sizeof(line)); //xlow xhi
        this->File->getline(line,sizeof(line)); //ylow yhi
        this->File->getline(line,sizeof(line)); //zlow zhi
        this->File->getline(line,sizeof(line)); //Item Entries
    }

    points = vtkPoints::New();
    points->SetDataTypeToFloat();
    points->Reset();

    double x1[3],x2[3];
    double F[3],N[3],Fn[3],Fs[3];
    double H[6],C[5];
    double N_mag,Fmag;
    int lc = 0;
    int pc = 0;

    std::string item;
    //printf("expecting %d elements\n",COUNT);

    // Allocate memory
    int *id1 = new int[COUNT];
    if (id1 == NULL)
    {
        printf("Error allocating memory for id1!\n");
        return 0;
    }
    int *id2 = new int[COUNT];
    if (id2 == NULL)
    {
        printf("Error allocating memory for id2!\n");
        return 0;
    }

    // Setup the force array
    vtkDoubleArray *forces = vtkDoubleArray::New();
    forces->SetNumberOfComponents(3);
    forces->SetName("F");

    vtkDoubleArray *normals = vtkDoubleArray::New();
    normals->SetNumberOfComponents(3);
    normals->SetName("N");

    vtkDoubleArray *shear = vtkDoubleArray::New();
    shear->SetNumberOfComponents(3);
    shear->SetName("shear");

    vtkDoubleArray *nforce = vtkDoubleArray::New();
    nforce->SetNumberOfComponents(3);
    nforce->SetName("nforce");

    vtkDoubleArray *hist = vtkDoubleArray::New();
    hist->SetNumberOfComponents(6);
    hist->SetName("history");

    vtkDoubleArray *cp = vtkDoubleArray::New();
    cp->SetNumberOfComponents(7);
    cp->SetName("ContactProperties");
    cp->SetComponentName(0,"Area");
    cp->SetComponentName(1,"ri");
    cp->SetComponentName(2,"rj");
    cp->SetComponentName(3,"distance");
    cp->SetComponentName(4,"overlap");
    cp->SetComponentName(5,"id1");
    cp->SetComponentName(6,"id2");

    vtkDoubleArray *rad = vtkDoubleArray::New();
    rad->SetNumberOfComponents(1);
    rad->Reset();
    rad->SetNumberOfTuples(2*COUNT);
    rad->SetComponentName(0,"radius");
    rad->SetName("Radius");

    vtkDoubleArray *vol = vtkDoubleArray::New();
    vol->SetNumberOfComponents(1);
    vol->Reset();
    vol->SetNumberOfTuples(COUNT);
    vol->SetComponentName(0,"Volume");
    vol->SetName("Void");

    vtkDoubleArray *pid = vtkDoubleArray::New();
    pid->SetNumberOfComponents(1);
    pid->Reset();
    pid->SetNumberOfTuples(2*COUNT);
    pid->SetComponentName(0,"ID");
    pid->SetName("PID");

    //Create a cell array to store the lines in and add the lines to it
    lines = vtkCellArray::New();
    H[0]=H[1]=H[2]=H[3]=H[4]=H[5]=0;
    F[0]=F[1]=F[2]=0;
    Fn[0]=Fn[1]=Fn[2]=Fs[0]=Fs[1]=Fs[2]=0;

    while ( this->File->getline(line,sizeof(line)) ) //every line
    {
        int ic=0;
        std::string line_content=line;
        std::stringstream ls(line_content);
        while(std::getline(ls, item, ' ')) //every item
        {
            switch (ic)
            {
            case 0:  x1[0]   = atof(item.c_str()); break;
            case 1:  x1[1]   = atof(item.c_str()); break;
            case 2:  x1[2]   = atof(item.c_str()); break;
            case 3:  x2[0]   = atof(item.c_str()); break;
            case 4:  x2[1]   = atof(item.c_str()); break;
            case 5:  x2[2]   = atof(item.c_str()); break;
            case 6:  id1[lc] = atoi(item.c_str()); break;
            case 7:  id2[lc] = atoi(item.c_str()); break;
            //8=periodic
            case 9:  F[0] = atof(item.c_str()); break;
            case 10: F[1] = atof(item.c_str()); break;
            case 11: F[2] = atof(item.c_str()); break;

            //5 contact-properties values
            case 12: C[0] = atof(item.c_str()); break; //contact Area
            case 13: C[1] = atof(item.c_str()); break; //ri
            case 14: C[2] = atof(item.c_str()); break; //rj
            case 15: C[3] = atof(item.c_str()); break; //d
            case 16: C[4] = atof(item.c_str()); break; //overlap distance

            //6 history values for hertz/history with rolling friction:
            case 17: H[0] = atof(item.c_str()); break;
            case 18: H[1] = atof(item.c_str()); break;
            case 19: H[2] = atof(item.c_str()); break;
            case 20: H[3] = atof(item.c_str()); break;
            case 21: H[4] = atof(item.c_str()); break;
            case 22: H[5] = atof(item.c_str()); break;
            }

            ++ic;
        } //item

        points->InsertNextPoint(x1[0], x1[1], x1[2]);
        points->InsertNextPoint(x2[0], x2[1], x2[2]);

        pid->InsertTuple1(pc,   id1[lc]);  // ID of point1
        pid->InsertTuple1(pc+1, id2[lc]);  // ID of point2

        C[5] = id1[lc];
        C[6] = id2[lc];

        rad->InsertTuple1(pc,   C[1]);  // radius of point1
        rad->InsertTuple1(pc+1, C[2]);  // radius of point2

        // see http://mathworld.wolfram.com/Sphere-SphereIntersection.html
        // V= PI*(Ri+Rj-d)^2*[d^2+2d*(Ri+Rj)-3*(Ri-Rj)^2]*/(12*d)
        double Ri = C[1];
        double Rj = C[2];
        double d  = C[3];
        double V  = M_PI*((Ri+Rj-d)*(Ri+Rj-d)*(d*d+2*d*Rj-3*Rj*Rj+2*d*Ri+6*Rj*Ri-3*Ri*Ri))/(12*d);

        N[0] = x2[0] - x1[0];
        N[1] = x2[1] - x1[1];
        N[2] = x2[2] - x1[2];

        N_mag = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
        N[0] = N[0] / N_mag;
        N[1] = N[1] / N_mag;
        N[2] = N[2] / N_mag;

        Fmag = sqrt(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);
        Fn[0] = Fmag*N[0];
        Fn[1] = Fmag*N[1];
        Fn[2] = Fmag*N[2];

        Fs[0] = F[0]-Fn[0];
        Fs[1] = F[1]-Fn[1];
        Fs[2] = F[2]-Fn[2];

        SHEAR += sqrt(Fs[0]*Fs[0] + Fs[1]*Fs[1] + Fs[2]*Fs[2]);

        tmpLine = vtkLine::New();
        tmpLine->GetPointIds()->SetId(0,pc);
        tmpLine->GetPointIds()->SetId(1,pc+1);
        lines->InsertNextCell(tmpLine);
        tmpLine->Delete();
        tmpLine = NULL;

        forces->InsertNextTupleValue(F);
        normals->InsertNextTupleValue(N);
        shear->InsertNextTupleValue(Fs);
        nforce->InsertNextTuple(Fn);
        hist->InsertNextTuple(H);
        cp->InsertNextTuple(C);
        vol->InsertTuple1(lc,V);

        pc += 2;
        ++lc;
    }//line

    // get the info object
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the ouptut
    vtkPolyData *myoutput = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    myoutput->SetPoints(points);
    myoutput->SetLines(lines);

    myoutput->GetPointData()->AddArray(pid);
    myoutput->GetPointData()->AddArray(rad);


    myoutput->GetCellData()->AddArray(forces);
    myoutput->GetCellData()->AddArray(normals);
    myoutput->GetCellData()->AddArray(shear);
    myoutput->GetCellData()->AddArray(nforce);
    myoutput->GetCellData()->AddArray(hist);
    myoutput->GetCellData()->AddArray(cp);
    myoutput->GetCellData()->AddArray(vol);

    // free memory
    points->Delete();
    lines->Delete();
    rad->Delete();
    vol->Delete();
    pid->Delete();

    forces->Delete(); forces = NULL;
    normals->Delete(); normals = NULL;
    shear->Delete(); shear = NULL;
    nforce->Delete(); nforce = NULL;
    hist->Delete(); hist = NULL;
    cp->Delete(); cp = NULL;


    vtkIntArray *intValue;
    intValue=vtkIntArray::New();
    intValue->SetNumberOfComponents(1);
    intValue->SetName("Dumpstep");
    intValue->InsertNextValue(TS);
    myoutput->GetFieldData()->AddArray(intValue);
    intValue->Delete();

    intValue = vtkIntArray::New();
    intValue->SetNumberOfComponents(1);
    intValue->SetName("COUNT");
    intValue->InsertNextValue(COUNT);
    myoutput->GetFieldData()->AddArray(intValue);
    intValue->Delete();

    vtkDoubleArray *doubleValue;
    doubleValue = vtkDoubleArray::New();
    doubleValue->SetNumberOfComponents(1);
    doubleValue->SetName("SHEARSUM");
    doubleValue->InsertNextValue(SHEAR);
    myoutput->GetFieldData()->AddArray(doubleValue);
    doubleValue->Delete();

    delete [] id1;
    delete [] id2;

    return 1;
}

int liggghts_forcechainreader::RequestInformation(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    //outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),-1);
    outInfo->Set(CAN_HANDLE_PIECE_REQUEST(),
                 1);
    return 1;
}


void liggghts_forcechainreader::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os  << indent << "FileName: "
        << (this->FileName ? this->FileName : "(NULL)") << endl;
}

int liggghts_forcechainreader::CanReadFile(const char *fname)
{
    return 1;
}
