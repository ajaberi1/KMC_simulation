
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>


using namespace std;

vector<int> find_nearest_neighbor(vector< vector<string> > atoms_list, vector<string> selected_atom, double niegbour_distance){
    // cout<<"selected atoms: "<<selected_atom[5]<<endl;
    vector<int> neighbors_index;
    double x = stod(selected_atom[1]);
    double y = stod(selected_atom[2]);
    double z = stod(selected_atom[3]);
    double index = stod(selected_atom[5]);
    for (int i = 0; i< atoms_list.size(); i++) {
        vector<string> loop_atom = atoms_list[i];
        string symbol_loop_atom = loop_atom[0];
        double x_loop_atom = stod(loop_atom[1]);
        double y_loop_atom = stod(loop_atom[2]);
        double z_loop_atom = stod(loop_atom[3]);
        int index_loop_atom = stod(loop_atom[5]);
        double distance;
        string Li = "Li";
        string X = "X";
        if (symbol_loop_atom == Li || symbol_loop_atom == X){
            distance = sqrt(pow((x-x_loop_atom),2) + pow((y-y_loop_atom),2) + pow((z-z_loop_atom),2));
            if (distance < niegbour_distance && distance != 0){
                neighbors_index.push_back(index_loop_atom);
            }
        }
    }
    return neighbors_index;
}


int main(){

    vector< vector<string> > atoms_list;
    ifstream myfile;
    string line;
    
    int num_row = 0;
    myfile.open("atoms_list_c_code.csv");
    while (getline(myfile, line)){
        num_row++;
    }
    myfile.close(); 

    vector<string> row;
    myfile.open("atoms_list_c_code.csv");
    string new_line;

    int count = 1;
    while (count <= num_row)
    {   
        count++;
        getline(myfile,new_line,'\n');
        stringstream sstr(new_line);
        vector<std::string> v;
        while(sstr.good())
        {
            string substr;
            getline(sstr, substr, ',');
            v.push_back(substr);
        }
        for (std::size_t i = 0; i < v.size(); i++){
            // cout << v[i] << endl;
            row.push_back(v[i]);
        }
        atoms_list.push_back(row);
        for (std::size_t i = 0; i < v.size(); i++){
            // cout << v[i] << endl;
            row.pop_back();
        }
            
    }
    myfile.close();

    for (int i = 0; i< atoms_list.size(); i++) {
        vector<string> atom = atoms_list[i];
        string symbol = atom[0];
        float x = stod(atom[1]);
        float y = stod(atom[2]);
        float z = stod(atom[3]);
        float index = stod(atom[5]);
        // cout<<symbol<<" "<<x<<" "<<y<<" "<<z<<" "<<index<<endl;
    }
  

    ofstream nieghbors_file;
    // NOTE: In this neighbors_c_code.csf file the selected atom in printed in the single row, and its neighbors are printed in the next row
    nieghbors_file.open ("neighbors_c_code.csv");
    for (int i = 0; i< atoms_list.size(); i++) {
        vector<string> selected_atom = atoms_list[i];
        string symbol = selected_atom[0];
        double x = stod(selected_atom[1]);
        double y = stod(selected_atom[2]);
        double z = stod(selected_atom[3]);
        string index_selected_atom = selected_atom[5];
        string Li = "Li";
        string X = "X";
        if (symbol == Li || symbol == X){
            vector<int> nearest_neighbor = find_nearest_neighbor(atoms_list, selected_atom, 3);

            // In this part we make a string to write in the .csv file. It has this format: index_selected_atom,nearest_neighbor[0],nearest_neighbor[1],....,nearest_neighbor[6]\n
            // This is an example for nieghbors_row: "2341,234,231,256\n" (the first one is the index of the selected atom)
            nieghbors_file << "Selected_atom:," + index_selected_atom + "\n";
            string nieghbors_row;
            nieghbors_row = "";
            for (int j = 0; j< nearest_neighbor.size(); j++) {
                nieghbors_row += to_string(nearest_neighbor[j]);
                nieghbors_row += ",";
            }
            if (nieghbors_row.size()>1){
                // Tp remove the extra Camma that we have from the previous For loop
                nieghbors_row.pop_back();
            }
            nieghbors_row+= "\n";
            // cout<<index_selected_atom;
            nieghbors_file << nieghbors_row;
            }
        }
    nieghbors_file.close();
}
