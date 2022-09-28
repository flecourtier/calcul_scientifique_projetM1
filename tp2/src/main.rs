mod facto;
use facto::*; 

use std::f64::consts::PI;

use std::fs::File;
use std::io::BufWriter;
use std::io::{Write};

//on définit u_0
fn u_0(x: f64) -> f64 {
    if x>=3./8. && x<=5./8. {
        return 1.;
    }
    else {
        return 0.;
    }
}

//on choisit L=1 pour les applications
//nous rentrons directement les valeurs sans utiliser la fonction u_0
fn coeff_ci(i: i32) -> f64 {
    if i==0 {
        return 1./4.;
    }
    else{
        let pi = PI;
        let sin = f64::sin;
        return 2./(pi*i as f64)*(sin(pi*5./8.*i as f64)-sin(pi*3./8.*i as f64));
    } 
}

//on définit les fonctions f et g
fn f_i(i: i32, x: f64)-> f64 {
    let pi = PI;
    let cos = f64::cos;
    return cos((i as f64)*pi*x);
}

fn g_i(i: i32, t: f64)-> f64 {
    let pi2 = PI*PI;
    let exp = f64::exp;
    return exp(-((i*i) as f64)*pi2*t);
}

// on calcule la somme
fn u_N(x: f64, t: f64, N: i32) -> f64 { //=somme c_i*g_i(t)*f_i(x)
    let mut somme=0.;
    for i in 0..(N+1) {
        somme=somme+coeff_ci(i)*g_i(i,t)*f_i(i,x);
    }
    return somme;
}

// affichage solution par méthode de développement en série
fn plotu_dvmpt_serie(xn: &Vec<f64>, tab_t: &Vec<f64>, N: i32) -> std::io::Result<()> {
    let fic = File::create("resu_serie_N500.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...
    
    let nb_x = xn.len();
    
    //pour écrire dans notre fichier la légende (nom des colonnes)
    write!(fic,"x ")?;
    for i in 0..tab_t.len() {
        write!(fic,"t={:?} ",tab_t[i])?;
    }
    writeln!(fic)?;

    for i in 0..nb_x {
        write!(fic,"{:?} ",xn[i])?;
        for t in tab_t {
            write!(fic,"{:?} ",u_N(xn[i],*t,N))?;
        }
        writeln!(fic)?;
    }

    Ok(())
}

//QUESTION 5 : theta schéma

//définit dt en fonction de theta
fn define_dt(t: f64, theta: f64, dx: f64) -> f64 {
    if theta>=1./2. {
        return t*dx;
    }
    else {
        return dx*dx/2./(1.-2.*theta);
    }
}

//résout l'équation
fn theta_schema(theta: f64, nb_x: usize, t_max: f64) -> (Vec<f64>, Vec<f64>) {
    let dx = 1./(nb_x as f64);
    let mut x_n = vec![0.;nb_x];
    for i in 0..nb_x {
        x_n[i] = (i as f64)*dx+dx/2.;
    }
    let dt=define_dt(t_max,theta,dx);

    // on définit la matrice M
    let low_M = vec![-1.*theta*dt/dx/dx;nb_x-1];
    let mut diag_M = vec![1.+2.*theta*dt/dx/dx; nb_x];
    diag_M[0] = 1.+1.*theta*dt/dx/dx;
    diag_M[nb_x-1] = 1.+1.*theta*dt/dx/dx;   
    let sup_M = vec![-1.*theta*dt/dx/dx;nb_x-1];

    // on définit la matrice N
    let low_N = vec![(1.-theta)*dt/dx/dx;nb_x-1];
    let mut diag_N = vec![1.-2.*(1.-theta)*dt/dx/dx; nb_x];
    diag_N[0] = 1.-(1.-theta)*dt/dx/dx;
    diag_N[nb_x-1] = 1.-(1.-theta)*dt/dx/dx;   
    let sup_N = vec![(1.-theta)*dt/dx/dx;nb_x-1];  

    //on détermine les factorisations LU de M et N
    let (low_M,diag_M,sup_M) = factolu(low_M,diag_M,sup_M);
    let (low_N,diag_N,sup_N) = factolu(low_N,diag_N,sup_N);
    //a noter que la fonction produit_matvec a besoin d'une matrice sous forme LU

    // on définit U_0
    let mut U_n = vec![0.;nb_x];
    for i in 0..nb_x {
        U_n[i]=u_0(x_n[i]);
    }

    // on résout le schéma en temps
    let mut t=0.;  
    let mut b = vec![0.;nb_x];
    while f64::abs((t_max as f64)-t)>1e-6 {
        // on calcule N*U_n
        b = produit_matvec(&low_N, &diag_N, &sup_N, &U_n, b);
        U_n = lusolve(&low_M, &diag_M, &sup_M, U_n, &b);
        t=t+dt;
    }

    (x_n,U_n)
}

fn plotu_theta_schema(theta: f64, nb_x: usize, tab_t: &Vec<f64>) -> std::io::Result<()> {
    let fic = File::create("resu_theta_schema_N500.dat")?;
    let mut fic = BufWriter::new(fic);
    
    //pour écrire dans notre fichier la légende (nom des colonnes)
    //et stocker les vecteurs solutions dans un vecteur 2D
    let mut Un : Vec<Vec<f64>> = vec![vec![0.;nb_x];tab_t.len()];
    write!(fic,"x ")?;
    for i in 0..tab_t.len() {
        write!(fic,"t={:?} ",tab_t[i])?;
        Un[i]=theta_schema(theta, nb_x, tab_t[i]).1;
    }
    writeln!(fic)?;
    
    let x_n=theta_schema(theta, nb_x, 0.).0;
    for i in 0..nb_x {
        write!(fic,"{:?} ",x_n[i])?;
        for t in 0..tab_t.len() {
            write!(fic,"{:?} ",Un[t][i])?;
        }
        writeln!(fic)?;
    }

    Ok(())
}

fn main() {
    // Resolution par méthode de développement en série
    let tab_t = vec![0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1.];
    let N=500;
    let nb_x = 100;
    let dx = 1./(nb_x as f64);
    
    let mut x_n = vec![0.;nb_x];
    for i in 0..nb_x {
        x_n[i] = (i as f64)*dx;
    }

    plotu_dvmpt_serie(&x_n, &tab_t, N);

    // Resolution par méthode des différences finies 
    // avec intégration en temps par un theta-schéma
    // let tab_t = vec![0.0001,0.001,0.01,0.1,1.];
    // let theta=0.;
    // let nb_x = 500;
    // let t_max=1.;

    // plotu_theta_schema(theta,nb_x,&tab_t);

}
