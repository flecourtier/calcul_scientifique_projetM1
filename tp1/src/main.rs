mod solve_laplace;
use solve_laplace::*;

use std::f64::consts::PI;

use std::fs::File;
use std::io::BufWriter;
use std::io::{Write};

fn factolu(mut l: Vec<f64>, mut d: Vec<f64>, u: Vec<f64>) 
-> (Vec<f64>, Vec<f64>, Vec<f64>){
    let n = d.len();
    for i in 0..n - 1 {
        l[i] = l[i] / d[i];
        d[i + 1] = d[i + 1] - l[i] * u[i];
    }
    (l,d,u)
}

fn lusolve(l: &Vec<f64>, d: &Vec<f64>, u: &Vec<f64>, mut x: Vec<f64>, b: &Vec<f64>)
-> Vec<f64> {
    x[0] = b[0];
    let n = d.len();
    // descente
    for i in 1..n {
        x[i] = b[i] - l[i - 1] * x[i - 1];
    }
    //remontée
    x[n - 1] = x[n - 1] / d[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (x[i] - u[i] * x[i + 1]) / d[i];
    }
    x
}

// f(x) = -u''(x) + alpha u(x)       avec alpha=1
fn f(x: f64) -> f64 {
    let pi = PI;
    let cos = f64::cos;
    let powi = f64::powi;
    let y = cos(pi * x)
        * (0.9e1 * powi(cos(pi * x), 2) * pi * pi + powi(cos(pi * x), 2) - 0.6e1 * pi * pi);
    y
}

// u(x)=cos^3(pi*x)
fn u(x: f64) -> f64 {
    let pi = PI;
    let cos = f64::cos;
    let powi = f64::powi;
    let y = powi(cos(pi * x), 3);
    y
}

//u_exact=[u(x_0),...,u(x_n)]
pub fn u_exact(xn: &Vec<f64>) -> Vec<f64> {
    let n = xn.len();
    let mut uex = vec![0.; n];
    for i in 0..n {
        uex[i]=u(xn[i]);
    }
    uex
}

// uapp approché de u_exact (u_i)
pub fn solve_pb(xn: &Vec<f64>, L: f64, alpha : f64) -> Vec<f64> {
    // A noter que xn contient les points 0 et L (et que ce sont des inconnus)
    let n = xn.len(); //n = N+2
    let dx = L/(n as f64);

    // on définit la matrice tridiagonale
    let low = vec![-1./dx/dx;n-1];
    let mut diag = vec![2./dx/dx+alpha; n]; // n = N+2 inconnus
    diag[0] = 1./dx/dx+alpha;
    diag[n-1] = 1./dx/dx+alpha;
    let sup = vec![-1./dx/dx;n-1];

    // on définit le second membre du système Ax=b
    // b_i=f(x_i)
    let mut b = vec![1.; n];
    for i in 0..n {
        b[i]=f(xn[i]);
    }

    // calcul de la factorisation LU de la matrice
    let (low,diag,sup) = factolu(low,diag,sup);

    // résolution du système linéaire
    let uapp = vec![0.; n];
    let uapp = lusolve(&low,&diag,&sup,uapp,&b);

    uapp
}

// création d'un fichier erreur.dat où se trouve 
// l'erreur maximale en norme absolue en fonction de dx (échelle logarithmique)
// on teste avec plusieurs p et plusieurs lambda
pub fn ploterr(L: f64, alpha : f64) -> std::io::Result<()> {
    let fic = File::create("erreur_pb2.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...
    
    write!(fic,"dx erreur\n")?;

    // pour faire tendre dx->0   
    for i in 0..50 {
        let n = (i+1)*10;
        let dx = L/((n+2) as f64);

        let mut xn = vec![0.; n+2];
        for i in 0..n+2 {
            xn[i] = dx * i as f64 + dx / 2.;
        }

        let uapp=solve_pb(&xn,L,alpha);
        let uex=u_exact(&xn);

        let mut erreur = 0.;
        for j in 0..n+2 {
            if f64::abs(uex[j]-uapp[j])>erreur {
                erreur = f64::abs(uex[j]-uapp[j]);
            }
        }
        
        write!(fic,"{} {}\n", f64::ln(dx), f64::ln(erreur))?;
    }    

    Ok(())
}

pub fn plotu(xn: Vec<f64>, L: f64, alpha: f64) -> std::io::Result<()> {
    let fic = File::create("resu_pb2.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...
 
    let n = xn.len();
    let uapp=solve_pb(&xn, L, alpha);
    let uex=u_exact(&xn);
    for i in 0..n {
        writeln!(fic, "{} {} {}", xn[i], uapp[i],uex[i])?;
    }

    Ok(())
}

pub fn test_pb() -> std::io::Result<()> {
    let n = 100;
    const L: f64 = 1.;
    let alpha = 1.;

    let dx = L / ((n+2) as f64);

    // on crée notre discrétisation pour un certain n
    let mut xn = vec![0.; n+2];
    for i in 0..n+2 {
        xn[i] = dx * i as f64 + dx / 2.;
    }

    // on résout le problème
    let uapp=solve_pb(&xn,L,alpha);
    let uex=u_exact(&xn);
    let mut erreur = 0.;
    for i in 0..n+2 {
        if f64::abs(uex[i]-uapp[i])>erreur {
            erreur = f64::abs(uex[i]-uapp[i]);
        }
    }
    println!("erreur pour n={:?} : {:?}",n,erreur);
    plotu(xn,L,alpha);

    // on vérifie la convergence quand dx->0
    ploterr(L,alpha);

    Ok(())
}

fn main() -> std::io::Result<()> {
    test_pb()
}

// fn main() -> std::io::Result<()> {
//     test_laplace()
// }