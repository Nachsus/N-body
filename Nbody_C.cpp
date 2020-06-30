/*
Nbody code!
Props to Peter for teaching me C++
*/

#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <math.h>
#include <chrono>
#include <thread>
using namespace std;

double Pi = 3.14159265359;
double G = 1E-8;

double t = 4*10000;
int tsteps = 1;
double tstep = t/tsteps;

double start_box_xmin = -100;
double start_box_xmax = 100;
double start_box_ymin = -100;
double start_box_ymax = 100;


class Body{
    public:
    double xpos;
    double ypos;
    double dx;
    double dy;
    double mass;
    double accx;
    double accy;
    double velx;
    double vely;
    double dvx;
    double dvy;
    double d;
    double soft;
    double accxbrute;
    double accybrute;


    Body(double x, double y, double m, double vx, double vy){
        xpos = x;
        ypos = y;
        velx = vx;
        vely = vy;
        mass = m;
        accx = 0;
        accy = 0;
    }


    void add_force_brute(Body& body){
        dx = xpos - body.xpos;
        dy = ypos - body.ypos;
        d = sqrt(dx*dx + dy*dy);
        if(d != 0){
            soft = 1;
            if(d<soft){
                d += soft;
                accxbrute += -(G * body.mass)/(d*d) * (dx/d);
                accybrute += -(G * body.mass)/(d*d) * (dy/d);
                cout << "hello" << endl;
            }
            else{
                accxbrute += -(G * body.mass)/(d*d) * (dx/d);
                accybrute += -(G * body.mass)/(d*d) * (dy/d);
            }
        }
    }


    void update_vel(){
        dvx = tstep * accx;
        dvy = tstep * accy;
        velx += dvx;
        vely += dvy;
    }


    void update_pos(){
        dx = tstep * velx;
        dy = tstep * vely;
        xpos += dx;
        ypos += dy;
    }
};


class Branch{
    public:
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    vector<Body> bodies_to_check;
    vector<Branch> sons;
    vector<Body> bodies;
    int N;
    double xCOM;
    double yCOM;
    vector<Body> bodies_in_area;
    double newx;
    double newy;
    double mass;
    double x;
    double y;
    double Mass;
    double box;
    double boy;
    double l;
    double d;
    double dx;
    double dy;
    double soft;


    Branch(double xmi, double xma, double ymi, double yma, vector<Body> btc){
        xmin = xmi;
        xmax = xma;
        ymin = ymi;
        ymax = yma;
        bodies_to_check = btc;
        mass = 0;
        bodies_in_corners();
        if(N > 1){
            makebranches();
        }
        COM_func();
    }


    bool between(double xbody, double ybody){
        if(xbody > xmin && xbody < xmax && ybody > ymin && ybody < ymax){
            return true;
        }
        else{
            return false;
        }
    }


    void bodies_in_corners(){
        for (int i = 0; i < (int) bodies_to_check.size(); i++)
        {
            Body b = bodies_to_check[i];
            if(between(b.xpos,b.ypos)){
                bodies_in_area.push_back(b);
            }
        }
        bodies = bodies_in_area;
        N = bodies_in_area.size();
        
    }


    void makebranches(){
        newx = (xmin+xmax)/2;
        newy = (ymin+ymax)/2;
        Branch son1(xmin,newx,ymin,newy,bodies);
        Branch son2(xmin,newx,newy,ymax,bodies);
        Branch son3(newx,xmax,ymin,newy,bodies);
        Branch son4(newx,xmax,newy,ymax,bodies);
        sons.push_back(son1);
        sons.push_back(son2);
        sons.push_back(son3);
        sons.push_back(son4);
    }


    void COM_func(){
        if(N==0){
            mass = 0;
        }
        else{
            x = 0;
            y = 0;
            Mass = 0;
            for (Body &b: bodies)
            {
                x += b.xpos*b.mass;
                y += b.ypos*b.mass;
                Mass += b.mass;
            }
            
            xCOM = x/Mass;
            yCOM = y/Mass;
            mass = Mass;
        }
    }


    double theta(Body& body){
        box = body.xpos;
        boy = body.ypos;
        if(mass != 0){
            l = abs(xmin-xmax);
            d = sqrt( (xCOM-box)*(xCOM-box) + (yCOM-boy)*(yCOM-boy) );
            if(d==0){
                return 0;
            }
            else{
                return l/d;
            }
        }
        else{
            return 0;
        }
    }


    void open(Body& body){
        sons[0].walk(body);
        sons[1].walk(body);
        sons[2].walk(body);
        sons[3].walk(body);
    }


    void walk(Body& body){
        if(theta(body) > 0.3 && N > 1){
            open(body);
        }
        else{
            add_force(body);
        }
    }


    void add_force(Body& body){
        if(N > 0){
            dx = body.xpos - xCOM;
            dy = body.ypos - yCOM;
            d = sqrt( dx*dx + dy*dy );
            if(d != 0){
                soft = 1;
                if(d < soft){
                    d += soft;
                    body.accx += -(G * mass)/(d*d) * (dx/d);
                    body.accy += -(G * mass)/(d*d) * (dy/d);
                }
                else{
                    body.accx += -(G * mass)/(d*d) * (dx/d);
                    body.accy += -(G * mass)/(d*d) * (dy/d);
                }
            }
        }
    }
};


double randDouble(double min, double max) {
    double r = (max-min)*double(rand())/(double(RAND_MAX)+1)+min;
    return r;
}


vector<Body> body_gen_random(int N, double mass, double xmin,double xmax,double ymin,double ymax){
    vector<Body> all_bodies;
    for (int i=0; i<N; i++){
        all_bodies.push_back(Body(randDouble(xmin,xmax),randDouble(ymin,ymax),mass,0,0));
    }
    return all_bodies;
}


int main(){
    int N = 200;
    srand((unsigned int)time(NULL));
    vector<Body> all_bodies;
    all_bodies = body_gen_random(N,100000000000,start_box_xmin +10,start_box_xmax -10,start_box_ymin +10,start_box_ymax -10);
    for (int i = 0; i < tsteps; i++)
    {

        double start = clock();
        double duration;
        Branch br(start_box_xmin,start_box_xmax,start_box_ymin,start_box_ymax,all_bodies);
        for(Body &b: all_bodies){
            b.accx = 0;
            b.accy = 0;
            br.walk(b);
            
            b.accxbrute = 0;
            b.accybrute = 0;
            for(Body &bod: all_bodies){
                b.add_force_brute(bod);
            }
            
            b.update_vel();
            b.update_pos();
        }
        double end = clock();
        duration = end-start;
        cout << duration/ ((double) CLOCKS_PER_SEC) << endl;
    }
    
    vector<double> forces;
    for(Body &body: all_bodies){
        double f1 = abs( (body.accxbrute + body.accybrute)/2 );
        double f2 = abs( (body.accx + body.accy)/2 );
        forces.push_back(  abs( (f2-f1)/f1 )  );
    }
    cout << "----------" << endl;
    
    for (int i = 0; i < all_bodies.size(); i++)
    {
       cout << forces[i] << endl;
       cout << ',' << endl;
    }
    cout << "----------" << endl;
    

    vector<double> forcestest;
    for (int i = 0; i < all_bodies.size(); i++)
    {
        cout << (all_bodies[i].accx + all_bodies[i].accy)/2 << endl;
        cout << ',' << endl;
        forcestest.push_back(abs((all_bodies[i].accx + all_bodies[i].accy)/2));
    }
    cout << "----------" << endl;
    
}