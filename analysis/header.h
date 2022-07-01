float* getContentFromTH1(const TH1D& h)
{
    int Nx = h.GetNbinsX();
    float* content = new float[3*Nx+1];
        // [0:Nx] -> bin content (len = Nx)
        // [Nx:2*Nx] -> bin error (len = Nx)
        // [2*Nx:3*Nx+1] -> bin edges (len = Nx+1)
    for (int i = 0 ; i < Nx ; i++)
    {
        content[i]        = h.GetBinContent(i+1);
        content[Nx+i]     = h.GetBinError(i+1);
        content[2*Nx+i]   = h.GetBinLowEdge(i+1);
    }
    content[3*Nx] = h.GetBinLowEdge(Nx+1);
    
    return content;
}

float* getContentFromTH2(const TH2D& h)
{
    int Nx = h.GetNbinsX();
    int Ny = h.GetNbinsY();
    float* content = new float[2*Nx*Ny+Nx+Ny+2];
        // [0:Nx*Ny] -> bin content (len = Nx*Ny) : rows = y values, columns = x values
        // [Nx*Ny:2*Nx*Ny] -> bin error (len = Nx*Ny)
        // [2*Nx*Ny:2*Nx*Ny+Nx+Ny+2] -> bin edges (len = Nx+Ny+2)
    for (int x = 0 ; x < Nx; x++)
    {
        for (int y = 0 ; y < Ny; y++)
        {
            content[y + Ny*x] = h.GetBinContent(x+1,y+1);
            content[Nx*Ny + y + Ny*x] = h.GetBinError(x+1,y+1);
            if (x == 0)
                content[2*Nx*Ny + Nx + 1 + y] = h.GetYaxis()->GetBinLowEdge(y+1);
        }
        if (x == 0)
            content[2*Nx*Ny + Nx + Ny + 1] =  h.GetYaxis()->GetBinLowEdge(Ny+1);
        content[2*Nx*Ny + x] = h.GetXaxis()->GetBinLowEdge(x+1);
    }
    content[2*Nx*Ny + Nx] = h.GetXaxis()->GetBinLowEdge(Nx+1);
    return content;
}

float* getContentFromTH3(const TH3D& h)
{
    int Nx = h.GetNbinsX();
    int Ny = h.GetNbinsY();
    int Nz = h.GetNbinsZ();
    float* content = new float[2*Nx*Ny*Nz+Nx+Ny+Nz+3];
        // [0:Nx*Ny*Nz] -> bin content (len = Nx*Ny*Nz)
        // [Nx*Ny*Nz:2*Nx*Ny*Nz] -> bin error (len = Nx*Ny*Nz)
        // [2*Nx*Ny*Nz:2*Nx*Ny*Nz+Nx+Ny+Nz+3] -> bin edges (len = Nx+Ny+Nz+3)
    for (int x = 0 ; x < Nx; x++)
    {
        for (int y = 0 ; y < Ny; y++)
        {
            for (int z = 0 ; z < Nz; z++)
            {
                content[z + Nz*y + Nz*Ny*x] = h.GetBinContent(x+1,y+1,z+1);
                content[Nx*Ny*Nz + z + Nz*y + Nz*Ny*x] = h.GetBinError(x+1,y+1,z+1);
                if (x == 0)
                    content[2*Nx*Ny*Nz + Nx + Ny + 2 + z] = h.GetZaxis()->GetBinLowEdge(z+1);
            }
            if (x == 0){
                content[2*Nx*Ny*Nz + Nx + 1 + y] =  h.GetYaxis()->GetBinLowEdge(y+1);
                content[2*Nx*Ny*Nz + Nx + Ny + Nz + 2] = h.GetZaxis()->GetBinLowEdge(Nz+1);
            }
        }
        if (x == 0)
            content[2*Nx*Ny*Nz + Nx + Ny + 1] =  h.GetYaxis()->GetBinLowEdge(Ny+1);

        content[2*Nx*Ny*Nz + x] = h.GetXaxis()->GetBinLowEdge(x+1);
    }
    content[2*Nx*Ny*Nz + Nx] = h.GetXaxis()->GetBinLowEdge(Nx+1);
    return content;
}



TH1F fillTH1(const float* edges, const float* values, const float * errors, int N, std::string name)
{
    TH1F h = TH1F(name.c_str(),name.c_str(),N,edges);
    for (int i = 0 ; i < N ; i++)
    {
        h.SetBinContent(i+1,values[i]);
        h.SetBinError(i+1,errors[i]);
    }    
    return h;
}
TH1D fillTH1(const double* edges, const double* values, const double * errors, int N, std::string name)
{
    TH1D h = TH1D(name.c_str(),name.c_str(),N,edges);
    for (int i = 0 ; i < N ; i++)
    {
        h.SetBinContent(i+1,values[i]);
        h.SetBinError(i+1,errors[i]);
    }    
    return h;
}


TH2F fillTH2(const float* xedges, const float* yedges, const float* values, const float* errors, int Nx, int Ny, std::string name)
{
    TH2F h = TH2F(name.c_str(),name.c_str(),Nx,xedges,Ny,yedges);
    for (int x = 0 ; x < Nx ; x++)
    {
        for (int y = 0 ; y < Ny ; y++)
        {
            h.SetBinContent(x+1,y+1,values[y+x*Ny]);
            h.SetBinError(x+1,y+1,errors[y+x*Ny]);
        }
    }    
    return h;
}

TH2D fillTH2(const double* xedges, const double* yedges, const double* values, const double* errors, int Nx, int Ny, std::string name)
{
    TH2D h = TH2D(name.c_str(),name.c_str(),Nx,xedges,Ny,yedges);
    for (int x = 0 ; x < Nx ; x++)
    {
        for (int y = 0 ; y < Ny ; y++)
        {
            h.SetBinContent(x+1,y+1,values[y+x*Ny]);
            h.SetBinError(x+1,y+1,errors[y+x*Ny]);
        }
    }    
    return h;
}

TH3F fillTH3(const float* xedges, const float* yedges, const float* zedges, const float* values, const float* errors, int Nx, int Ny, int Nz, std::string name)
{
    TH3F h = TH3F(name.c_str(),name.c_str(),Nx,xedges,Ny,yedges,Nz,zedges);
    for (int x = 0 ; x < Nx ; x++)
    {
        for (int y = 0 ; y < Ny ; y++)
        {
            for (int z = 0 ; z < Nz ; z++)
            {
                h.SetBinContent(x+1,y+1,z+1,values[z + Nz*y + Nz*Ny*x]);
                h.SetBinError(x+1,y+1,z+1,errors[z + Nz*y + Nz*Ny*x]);
            }
        }
    }    
    return h;
}

TH3D fillTH3(const double* xedges, const double* yedges, const double* zedges, const double* values, const double* errors, int Nx, int Ny, int Nz, std::string name)
{
    TH3D h = TH3D(name.c_str(),name.c_str(),Nx,xedges,Ny,yedges,Nz,zedges);
    for (int x = 0 ; x < Nx ; x++)
    {
        for (int y = 0 ; y < Ny ; y++)
        {
            for (int z = 0 ; z < Nz ; z++)
            {
                h.SetBinContent(x+1,y+1,z+1,values[z + Nz*y + Nz*Ny*x]);
                h.SetBinError(x+1,y+1,z+1,errors[z + Nz*y + Nz*Ny*x]);
            }
        }
    }    
    return h;
}
