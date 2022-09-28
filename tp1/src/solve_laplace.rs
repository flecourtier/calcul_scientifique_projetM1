mod facto;
use facto::*; 

use std::f64::consts::PI;
use std::f64;
use std::fs::File;
use std::io::BufWriter;
use std::io::{Write};

/*
Applications : problème de Laplace en dimension 1
*/

//1er problème (Dirichlet)

//f(x)=-u''(x)
pub fn f(x : f64, lambda : f64, p : i32) -> f64 {
    let pi = PI;
    let cos = f64::cos;
    let sin = f64::sin;
    let powi = f64::powi;
    let exp = f64::exp;
    let p:f64 = p as f64;
    return exp(-lambda*x)*((powi(pi*p,2)-powi(lambda,2))*sin(pi*p*x)+2.*lambda*pi*p*cos(pi*p*x));
}

//u(x)=sin(pi*p*x) * exp(-lambda*x)
pub fn u(x: f64, lambda : f64, p : i32) -> f64 {
    let pi = PI;
    let sin = f64::sin;
    let exp = f64::exp;
    let p:f64 = p as f64;
    return sin(pi*p*x) * exp(-lambda*x);
}

//u_exact=[u(x_0),...,u(x_n)]
pub fn u_exact(xn: &Vec<f64>, lambda : f64, p : i32) -> Vec<f64> {
    let n = xn.len();
    let mut uex = vec![0.; n];
    for i in 0..n {
        uex[i]=u(xn[i],lambda,p);
    }
    uex
}

// uapp approché de u_exact (u_i)
pub fn solve_laplace(xn: &Vec<f64>, L: f64, lambda : f64, p : i32) -> Vec<f64> {
    // A noter que xn contient les points 0 et L (mais ne sont pas inconnus)
    let n = xn.len(); //n = N+2
    let dx = L/((n-1) as f64);

    // on définit la matrice tridiagonale
    let low = vec![-1./dx/dx;n-3];
    let mut diag = vec![2./dx/dx; n-2]; // N inconnus
    let sup = vec![-1./dx/dx;n-3];

    // on définit le second membre du système Ax=b
    // b_i=f(x_i)
    let mut b = vec![1.; n-2];
    for i in 0..n-2 {
        b[i]=f(xn[i+1],lambda,p);
    }

    // calcul de la factorisation LU de la matrice
    let (low,diag,sup) = factolu(low,diag,sup);

    // résolution du système linéaire
    let uapp = vec![0.; n-2];
    let uapp = lusolve(&low,&diag,&sup,uapp,&b);

    //on rajoute les conditions aux bords (N+2 points)
    let mut uapp_cond = vec![0.; n];
    for i in 1..n-1 {
        uapp_cond[i] = uapp[i-1];
    }

    uapp_cond
}

// création d'un fichier erreur.dat où se trouve 
// l'erreur maximale en norme absolue en fonction de dx (échelle logarithmique)
// on teste avec plusieurs p et plusieurs lambda
pub fn ploterr(L: f64, tab_lambda : &Vec<f64>, tab_p : &Vec<i32>) -> std::io::Result<()> {
    let fic = File::create("erreur_pb1.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...
    
    write!(fic,"dx ")?;

    //pour écrire dans notre fichier la légende (nom des colonnes)
    for l in 0..tab_lambda.len() {
        for p in 0..tab_p.len() {
            write!(fic,"lambda={:?},p={:?} ",tab_lambda[l],tab_p[p])?;
        }
    }
    writeln!(fic)?;

    // pour faire tendre dx->0   
    for i in 0..50 {
        let n = (i+1)*10;
        let dx = L/((n+1) as f64);

        let mut xn = vec![0.; n+2];
        for i in 0..n+2 {
            xn[i] = dx * (i as f64);
        }
        
        write!(fic, "{} ", f64::ln(dx))?;

        for l in 0..tab_lambda.len() {
            for p in 0..tab_p.len() {   

                let uapp=solve_laplace(&xn,L,tab_lambda[l],tab_p[p]);
                let uex=u_exact(&xn,tab_lambda[l],tab_p[p],);
        
                let mut erreur = 0.;
                for j in 0..n+2 {
                    if f64::abs(uex[j]-uapp[j])>erreur {
                        erreur = f64::abs(uex[j]-uapp[j]);
                    }
                }
                
                write!(fic,"{} ", f64::ln(erreur))?;

            }
        }

        writeln!(fic)?;
    }    

    Ok(())
}

pub fn plotu(xn: Vec<f64>, L: f64, lambda : f64, p : i32) -> std::io::Result<()> {
    let fic = File::create("resu_pb1.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...
 
    let n = xn.len();
    let uapp=solve_laplace(&xn,L,lambda,p);
    let uex=u_exact(&xn,lambda,p);
    for i in 0..n {
        writeln!(fic, "{} {} {}", xn[i], uapp[i],uex[i])?;
    }

    Ok(())
}

pub fn test_laplace() -> std::io::Result<()> {
    let n = 100;
    const L: f64 = 1.;
    let lambda=2.5;
    let p=4;

    let dx = L/((n+1) as f64);

    // on crée notre discrétisation pour un certain n
    let mut xn = vec![0.; n+2];
    for i in 0..n+2 {
        xn[i] = dx * (i as f64);
    }

    // on résout le problème
    let uapp=solve_laplace(&xn,L,lambda,p);
    let uex=u_exact(&xn,lambda,p);
    let mut erreur = 0.;
    for i in 0..n+2 {
        if f64::abs(uex[i]-uapp[i])>erreur {
            erreur = f64::abs(uex[i]-uapp[i]);
        }
    }
    println!("erreur pour n={:?} : {:?}",n,erreur);
    plotu(xn,L,lambda,p);

    // on vérifie la convergence quand dx->0
    let tab_lambda=vec![1.,2.,3.];
    let tab_p=vec![1,2,3];
    ploterr(L, &tab_lambda, &tab_p);

    Ok(())
}