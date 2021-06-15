#include <iostream>
#include <vector>
#include <iomanip>
#include <math.h>
#define db float
using namespace std;
extern "C" {void findpores(
    int, \
    int, \
    int&, db[], db[], db[], int[], \
    int&, db[], db[], db[],\
    int,  db[], db[], db[],\
    db, db, db, db, db, \
    db[], db[], db[]);}
//std::vector<db> linspace(db a, db b, int n);

inline std::vector<db> linspace(db a, db b, int n) {
    std::vector<db> array;
    db step = (b-a) / (n-1);
    for(int i=0; i<n; i++) 
        array.push_back(a + step * i); 
    return array;
}

inline db cal_min_d(db x1,db x2, db x3, int natom, db *xred1, db *xred2, db *xred3,db r00,db r01, db r02, db r10, db r11, db r12, db r20,db r21, db r22   ) {
    db dsqrmin = 100, d1, d2, d3, dx, dy, dz, dsqr ;
    //cout <<"natom "<<natom<<endl;
    for(int i_at=0; i_at<natom; i_at++) {
        //cout << xred1[n]<<endl;
        //Periodic boundary conditions
        d1 = x1 - xred1[i_at];
        if (d1 >  0.5) d1 = d1 - 1;
        if (d1 < -0.5) d1 = d1 + 1;


        d2 = x2 - xred2[i_at];
        if (d2 >  0.5) d2 = d2 - 1;
        if (d2 < -0.5) d2 = d2 + 1;



        d3 = x3 - xred3[i_at];
        if (d3 >  0.5) d3 = d3 - 1;
        if (d3 < -0.5) d3 = d3 + 1;




        dx = r00 * d1 + r10 * d2 + r20 * d3;
        dy = r01 * d1 + r11 * d2 + r21 * d3;
        dz = r02 * d1 + r12 * d2 + r22 * d3;
        dsqr = dx * dx + dy * dy + dz * dz;
        if (dsqr < dsqrmin) dsqrmin = dsqr; //find distance between current point and nearest atom
        //cout << dsqr<<" ";
        //cout <<"i_at "<<i_at<<endl;
        //cout <<xred1[i_at]<<" "<<xred2[i_at]<<" "<<xred3[i_at]<< endl;


    }
    return dsqrmin;
}





void findpores( int check_pore_vol, \
                int max_npores, \
                int &ntot,  db l_pxred1[], db l_pxred2[], db l_pxred3[], int l_npores[], \
                int &npores, db pxred1[],   db pxred2[],   db pxred3[],   \
                int natom,   db ixred1[],   db ixred2[],   db ixred3[],\
                db r_matrix, db r_impurity, db step_dec,   db fine, db prec, \
                db r1[3], db r2[3], db r3[3]) {
    //n_tot - total number of local points inside pores; the dimension of pxred arrays
    //npores - total number of pores added in xred after natom
    //db fine//Increase mesh  during local scanning

                    //Find all nearest points up to the surfaces of nearest spherical atoms.
                    //To prevent accounting of points we control rate of growth
                    //We assume that if the rate of growth starts to increase after decreasing, this means coalescence with neihgbour pore

    //
    db dsqrmin,prate,rate; 
    int iterat, i_porepoints, n_porepoints = 0, n_localpoints, natom_init = natom, i_pores = 0, i_tot = 0;
    bool reincreasing, decreasing; 
	std::vector<db> steps1, localsteps1;
	std::vector<db> steps2, localsteps2;
	std::vector<db> steps3, localsteps3;

    db r00 = r1[0], r01 = r1[1], r02 = r1[2];
    db r10 = r2[0], r11 = r2[1], r12 = r2[2];
    db r20 = r3[0], r21 = r3[1], r22 = r3[2];
    int nsteps1 = ceil( sqrt(r00*r00 + r01*r01 +r02*r02) / step_dec); // number of steps in each direction in reduced space
    int nsteps2 = ceil( sqrt(r10*r10 + r11*r11 +r12*r12) / step_dec);
    int nsteps3 = ceil( sqrt(r20*r20 + r21*r21 +r22*r22) / step_dec);
    db scans1f = 1. / nsteps1 / fine; db scans2f = 1. / nsteps2 / fine; db scans3f = 1. / nsteps3 / fine; //scanning fine step
   //cout <<"scansif "<<scansif<<endl;
    cout <<  "Number of main points = "<< nsteps1*nsteps2*nsteps3 <<endl;
	steps1 = linspace(0, 1, nsteps1); 
	steps2 = linspace(0, 1, nsteps2);
	steps3 = linspace(0, 1, nsteps3);

    for(int i=0; i<natom; i++) { //arrays for internal use, including coordinates of pores 
        pxred1[i] = ixred1[i];
        pxred2[i] = ixred2[i];
        pxred3[i] = ixred3[i];
    }

    db rmsqr  = r_matrix*r_matrix + step_dec*step_dec; // for small local balls to fill the pore
    db refsqr = (r_matrix + r_impurity)*(r_matrix + r_impurity); // square of distance between atoms and pores
    cout << std::setprecision(3);
    //cout <<xred1[n]<< endl;
    //return;
    db x1s, x2s, x3s, x1a, x2a, x3a, diff;
    prec = prec*prec;
    int npoints;

    #pragma omp parallel
    {
    #pragma omp for schedule(dynamic, 100)
    for( auto &x1 : steps1 ) 
        for( auto &x2 : steps2 )
            for( auto &x3 : steps3 ) {
                dsqrmin = cal_min_d(x1, x2, x3, natom, pxred1, pxred2, pxred3, r00,r01,r02,r10,r11,r12,r20,r21,r22 );
                //cout << "dsqrmin C "<<dsqrmin << endl;
                if (dsqrmin > refsqr) { //this point situated in the spherical pore with radius larger than ri, considering spherical matrix atoms 
                    //cout <<  "Pore was found\n";

                    //cout << "dsqrmin main "<<dsqrmin << endl;
                    //cout <<  "i, j, k "<<i<<" "<<j<<" "<<k<<endl;

                    
                    cout <<  "\n\nStart determine center\n  ";
                    iterat = int(10*fine);
                    npoints = 1;
                    //Determine precisely center of pore
                    while (1) { //Start determine center
                        x1a = x1; // save current center
                        x2a = x2; 
                        x3a = x3;
                        cout <<  "x1, x2, x3  "<< x1 <<' '<< x2 <<' ' << x3 <<endl;
                        x1s = x1; //sum
                        x2s = x2; 
                        x3s = x3;
                        npoints = 1; //number of points to average

                        localsteps1 = linspace( -iterat * scans1f, iterat * scans1f, iterat*2);
                        localsteps2 = linspace( -iterat * scans2f, iterat * scans2f, iterat*2);
                        localsteps3 = linspace( -iterat * scans3f, iterat * scans3f, iterat*2);

                        for( auto &lx1 : localsteps1 ) 
                            for( auto &lx2 : localsteps2 )
                                for( auto &lx3 : localsteps3 ) {
                                    //if (lx1*lx1 + lx2*lx2 +  lx3*lx3 > refsqr)
                                    dsqrmin = cal_min_d(x1+lx1, x2+lx2, x3+lx3, natom, pxred1, pxred2, pxred3,r00,r01,r02,r10,r11,r12,r20,r21,r22 );

                                    if (dsqrmin > refsqr) {
                                        x1s = x1s + (x1 + lx1);
                                        x2s = x2s + (x2 + lx2);
                                        x3s = x3s + (x3 + lx3);

                                        npoints++;
                                        
                                    }
                                }
                        cout <<  "number of points to average  "<< npoints <<endl;
                        x1 = x1s/npoints; //New averaged center
                        x2 = x2s/npoints;
                        x3 = x3s/npoints;
                        diff  = (x1-x1a)*(x1-x1a) + (x2-x2a)*(x2-x2a) + (x3-x3a)*(x3-x3a); //difference between previus and current center.
                        cout <<  "x1, x2, x3  "<< x1 <<' '<< x2 <<' ' << x3 <<endl;
                        cout << "diff " << sqrt(diff);
                        if ( diff < prec ) break; //the precicion of center
                        iterat += 1;
                        cout <<  "\niterat  "<< iterat <<endl;
                    }//end determine center                    


                    iterat = 1; reincreasing = false; decreasing = false; prate = 0; rate = 0;  //reset state
                    while (check_pore_vol) { //Start checking volume
                        localsteps1 = linspace( -iterat * scans1f, iterat * scans1f, iterat*2);
                        localsteps2 = linspace( -iterat * scans2f, iterat * scans2f, iterat*2);
                        localsteps3 = linspace( -iterat * scans3f, iterat * scans3f, iterat*2);
                        n_localpoints = 8 * iterat * iterat * iterat; //Total number of local points
                        //n_localpoints = localstepsi.size();
                        //cout <<  "iterat * 2  "<< iterat * 2 <<endl;
                        i_porepoints = 0; //number of points inside the current pore
                        cout <<  "OK!!!\n";
                        for( auto &lx1 : localsteps1 ) 
                            for( auto &lx2 : localsteps2 )
                                for( auto &lx3 : localsteps3 ) {
                                    dsqrmin = cal_min_d(x1+lx1, x2+lx2, x3+lx3, natom, pxred1, pxred2, pxred3,r00,r01,r02,r10,r11,r12,r20,r21,r22 );
                                    //cout << "dsqrmin local "<<dsqrmin;
                                    if (dsqrmin > rmsqr) {
                                        //if ()
                                        //cout <<"N "<< n_porepoints<< " included "<<dsqrmin << endl;
                                        //cout <<  "li, lj, lk "<<li<<" "<<lj<<" "<<lk<<endl;
                                        l_pxred1[i_tot + i_porepoints] = x1 + lx1; //add local pores
                                        l_pxred2[i_tot + i_porepoints] = x2 + lx2;
                                        l_pxred3[i_tot + i_porepoints] = x3 + lx3;
                                        i_porepoints += 1;
                                        
                                    }
                                }

                        n_porepoints = i_porepoints;
                        prate = rate;
                        rate = (db) n_porepoints / n_localpoints;// # part from local points inside the pore for this iteration

                        //cout << "Number of local points "       << n_localpoints    << endl;
                        //cout << "Number of points inside pore " << n_porepoints     << endl;
                        //cout << "rate = "                       << rate             << endl;
                        //if (rate < prate-0.1) decreasing = true; //becomes true in the case of decrasing with high speed 
                        //if (decreasing and rate > prate-0.1) break;
                        //if (decreasing and fabs(rate - prate)<0.1 ) break;
                        //if (iterat > 5) break;
                        if ( rate <0.6 ) break;
                        iterat += 1;
            
                    }//end checking volume

                    pxred1[natom] = x1; //add coordinates of pore to arrays (it is equal adding matrix atom! to pore and not impurity atom)
                    pxred2[natom] = x2;
                    pxred3[natom] = x3;
                    natom++;
                    if (natom > max_npores) cout <<  "Error natom > than maxat\n";
                    i_tot += n_porepoints; //total number of points inside all pores
                    //Volume of pore:
                    //db a = step_dec/fine; //the side of little cube formed by the mesh which is used to find spheres inside the pore.
                    //cout <<  "Volume of pore is "<< n_porepoints * a*a*a <<" A^3\n";
                    if (i_tot > max_npores) {cout <<  "Error n_tot > than maxat. n_tot set to zero!\n"; i_tot = 0;}

                    l_npores[i_pores] = n_porepoints; //contains number of small balls in each pore.
                    cout <<"\nImpurity added; number = "         << i_pores       << endl;
                    i_pores+=1; 
                }
            }

    }


    ntot = i_tot;
    npores = i_pores;
    //cout <<"Initial number of atoms = " << natom_init   << endl;
    cout <<"\nNumber of pores = "         << npores       << endl;

    //for(int n=0; n<natom; n++) {
    //    cout << pxred1[n]<<" "<< pxred2[n]<<" "<< pxred3[n]<<endl;
    //}
    return;
}





