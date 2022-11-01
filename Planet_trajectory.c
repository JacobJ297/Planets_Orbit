/*
* Project (2) for the PHY2027 module
* Author: Jacob Johnson
* Date: 09/12/2021
*/

/*
 *
 * Program description: This program integrates the equations of planetary motion by reading in the planets initial conditions and then
 *                      numerically integrating with the second order Runga-Kutta method which is a time stepping method. It plots a graph
 *                      of the chosen planet's orbit for the user in gnuplot.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NAMELENGTH 20
#define LINELENGTH 250


typedef struct initial {
char name[NAMELENGTH];
float x0;
float y0;
float vx0;
float vy0;
}initial;

typedef struct constants{
float g;
float sunm;
float dt;
float n;
} constants;

void read(int planetnum, initial *p );
void rk2(initial *p, constants c);

int main()
{
    //Set up structures
    constants c={6.67E-11,1.989E30,86400,0}; //Will change dt and n depending on planet due to different times taken to orbit sun
    initial *pinitial=NULL;
    pinitial=(initial*)malloc(sizeof(initial));

    // menu to choose which planet to visualise
    int option;
    while(option != 9){
    do {
        printf("Menu- please choose a planets orbit to calculate and graph\n"
        "1. Mercury \n"
        "2. Venus \n"
        "3. Earth\n"
        "4. Mars\n"
        "5. Jupiter\n"
        "6. Saturn\n"
        "7. Uranus \n"
        "8. Neptune\n"
        "9. quit\n");
        scanf("%d", &option);

        if (option!=1 && option!=2 &&option!=3 &&option!=4 &&option!=5 &&option!=6 &&option!=7 &&option!=8 &&option!=9){
            printf("This option does not exist, please enter a correct option from the menu \n");
        }
    } while (option!=1 && option!=2 && option!=3 && option!=4 &&option!=5 &&option!=6 &&option!=7 &&option!=8 &&option!=9);

    switch (option) {


    case 1:
    printf("You should have a window of a graph of Mercury's orbit. Mercury takes 88 days to complete this orbit around the sun. \n");
    c.n=88;
    read(1,pinitial);
    rk2(pinitial,c);
    break;

    case 2:
    printf("You should have a window of a graph of Venus's orbit. Venus takes 225 days to complete this orbit around the sun. \n");
    c.n=250;
    read(2,pinitial);
    rk2(pinitial,c);
    break;

    case 3:
    printf("You should have a window of a graph of Earth's orbit. Earth takes 365.25 days to complete this orbit around the sun. \n");
    c.n=400;
    read(3,pinitial);
    rk2(pinitial,c);
    break;

    case 4:
    printf("You should have a window of a graph of Mars's orbit. Mars takes 687 days to complete this orbit around the sun. \n");
    c.n=700;
    read(4,pinitial);
    rk2(pinitial,c);
    break;

    case 5:
    printf("You should have a window of a graph of Jupiter's orbit. Jupiter takes 12 years to complete this orbit around the sun. \n");
    c.n=4500;
    read(5,pinitial);
    rk2(pinitial,c);
    break;

    case 6:
    printf("You should have a window of a graph of Saturn's orbit. Saturn takes 29 years to complete this orbit around the sun. \n");
    c.n=11000;
    read(6,pinitial);
    rk2(pinitial,c);
    break;

    case 7:
    printf("You should have a window of a graph of Uranus's orbit. Uranus takes 84 years to complete this orbit around the sun. \n");
    c.dt=1.6E7;
    c.n=170;
    read(7,pinitial);
    rk2(pinitial,c);
    break;

    case 8:
    printf("You should have a window of a graph of Neptune's orbit. Neptune takes 165 years to complete this orbit around the sun. \n");
    c.dt=1.6E7;
    c.n=350;
    read(8,pinitial);
    rk2(pinitial,c);
    break;

    case 9:
    printf("Goodbye, Thanks for using the menu");
    return 0;
    break;
           }
        }
    return 0;
}

//Function to read in initial conditions for planet
void read(int planetnum, initial *p ){
    char line[LINELENGTH];
    //Open file of initial conditions
    FILE *infile;
    infile = fopen("Initial_Conditions.txt", "r");
    //Discard top line and other planets before
    for (int i=0; i<planetnum; ++i){
    fscanf(infile, "%[^\n]\n", line);
    }
    //Read in the initial conditions and save to structure
    fscanf(infile, "%[^\n]\n", line);
    sscanf(line," %s %g %g %g %g", &p->name, &p->x0, &p->y0, &p->vx0, &p->vy0);
    return;
}

// Function to calculate planets orbit using rk2 and plot using gnuplot
 void rk2(initial *p,constants c){
    //Open file we will write results to
    FILE *outfile=NULL;
    outfile = fopen("orbit1.tmp", "w");

    //Setting up gnuplot for plotting graph
    FILE *gnupipe=NULL;
    gnupipe=_popen("gnuplot -persistent","w");
    // I will send these commands to gnuplot to plot graph
    char *Gnucommands [] ={"set style data lines","set title \"orbit\"", "plot 'orbit1.tmp'"};

    //Rename structure variables for easier use in calculations
    float x, y, vx, vy;
    x=p->x0;
    y=p->y0;
    vx=p->vx0;
    vy=p->vy0;
    double k21, k22, k23, k24, base, base2;

    // calculations
    for (int i=0; i<c.n; ++i){

        fprintf(outfile," %f \t %f \n", x, y);
        base=pow(x,2)+pow(y,2);
        base2=pow(x+0.5*c.dt*vx,2)+pow(y+0.5*c.dt*vy,2);
        k21=c.dt*(vx-(0.5*c.dt*c.g*c.sunm*x)/pow(base,1.5));
        k22=c.dt*(vy-(0.5*c.dt*c.g*c.sunm*y)/pow(base,1.5));
        k23=c.dt*((-c.sunm*c.g*(x+0.5*c.dt*vx))/(pow(base2,1.5)));
        k24=c.dt*((-c.sunm*c.g*(y+0.5*c.dt*vy))/(pow(base2,1.5)));

        x+=k21;
        y+=k22;
        vx+=k23;
        vy+=k24;
    }
    // plot to gnuplot
    for (int i =0; i<3; i++){
        fprintf(gnupipe,"%s\n", Gnucommands[i]);

    }

    fclose(outfile);
    fclose(gnupipe);
    return;
 }

/* Results from program. Please note this is for one planet mercury but each planet can be chosen through menu. The graph of Mercury's orbit is provided in zip file.

Menu- please choose a planets orbit to calculate and graph
1. Mercury
2. Venus
3. Earth
4. Mars
5. Jupiter
6. Saturn
7. Uranus
8. Neptune
9. quit
1
You should have a window of a graph of Mercury's orbit. Mercury takes 88 days to complete this orbit around the sun.
*/
