// VolEsti

// Copyright (c) 2012-2017 Vissarion Fisikopoulos

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

int linear_extensions_to_order_polytope(std::istream &is,
                                        std::ostream &os){

    std::string ns;
    //read n: the number of elements in the poset
    std::getline(is, ns, ' ');
    std::istringstream buffer(ns);
    int n;
    buffer >> n;
    //read m: the number of edges/relations in the poset
    std::getline(is, ns, '\n');
    std::istringstream bufferm(ns);
    int m;
    bufferm >> m;
    
    os << "order_"<<n<<".ine\n";
    os << "H-representation\n";
    os << "begin\n";
    os << " " << 2*n+m << " " << n+1 << " integer\n";

    //this is a cube
    for(int i=0; i<n; ++i){
        os << " " << 0 << " ";
        for(int j=0; j<n; ++j){
            if(i==j)
                os << 1 << " ";
            else os << 0 << " ";
        }
        os << "\n";
    }
    for(int i=0; i<n; ++i){
        os << " " << 1 << " ";
        for(int j=0; j<n; ++j){
            if(i==j)
                os << -1 << " ";
            else os << 0 << " ";
        }
        os << "\n";
    }
    
    //add constraints for the poset relations
    std::string point;
    while(!std::getline(is, point, ']').eof()) {
        std::vector<int> ipoint;
        point.erase( std::remove( point.begin(), point.end(), ' ' ),
                     point.end() );
        point.erase( std::remove( point.begin(), point.end(), '[' ),
                     point.end() );
        std::string coord;
        if (point[0]==',')
            point.erase(0,1);
        std::stringstream stream(point);
        if (!point.empty()){
            int i=1;
            os<<" 0 ";
            //read point 1
            std::getline(stream, coord, ',');
            std::istringstream buffer(coord);
            int temp;
            buffer >> temp;
            ipoint.push_back(temp);
            //os << temp << " ";
            //write inequality 1
            for(;i<temp;++i)
                os<<"0 ";
            os<<"1 ";++i;
            //read point 2
            std::getline(stream, coord, ',');
            std::istringstream buffer2(coord);
            int temp2;
            buffer2 >> temp2;
            ipoint.push_back(temp2);
            //os << temp2 << "\n";
            //write inequality 2
            for(;i<temp2;++i)
                os<<"0 ";
            os<<"-1 ";++i;
            //fill the rest with 0s
            for(;i<=n;++i)
                os<<"0 ";
            os<<"\n";
            //pointset.push_back(ipoint);
        }
    }
    os << "end\ninput_incidence" << std::endl;
}
