use std::f64::consts::PI;

use std::fs::File;
use std::io::BufWriter;
use std::io::{Write};

//QUESTION 3 : TRAFIC AUTOMOBILE

//u_0
fn u_0(x: f64, L: f64, pho0: u32) -> f64 {
    if (x>=0.) & (x<=L) {
        return pho0 as f64;
    }
    else {
        return 0.;
    }
}

//f(rho)=rho*v0*(1-pho/rho0)
fn f(u: f64, v0: u32, pho0: u32) -> f64 {
    u*(v0 as f64)*(1.-u/(pho0 as f64))
}

//f'(rho)=v0-2 v0/rho0*rho
fn f_prime(u: f64, v0: u32, pho0: u32) -> f64 {
    (v0 as f64) - 2.*(v0 as f64)/(pho0 as f64)*u
}

//retourne lambda 
fn lambda(uL: f64, uR: f64, v0: u32, pho0: u32) -> f64 {
    f64::max( f64::abs(f_prime(uL, v0, pho0)) , f64::abs(f_prime(uR, v0, pho0)) )
}

//retourne la valeur du flux
fn flux(uL: f64, uR: f64, v0: u32, pho0: u32) -> f64 {
    let lambda = lambda(uL, uR, v0, pho0);
    (f(uL, v0, pho0)+f(uR, v0, pho0))/2.-lambda/2.*(uR-uL)
}

fn solve_rusanov_trafic(N: usize, L: f64, T: f64, v0: u32, pho0: u32) -> (Vec<f64>,Vec<f64>,f64) {
    let dx = L/(N as f64);
    let mut t=0.;

    //on définit les points de discrétisation
    let mut x=vec![0.;N+2];
    for i in 0..N+2 {
        x[i]=(i as f64 - 1./2.)*dx;
    }

    // on définit u au temps 0
    let mut uapp=vec![0.;N+2];
    for i in 0..N+2 {
        uapp[i]=u_0(x[i],L,pho0);
    }

    // on définit un vecteur qui sera u au temps n+1
    let mut uapp_p1=vec![0.;N+2];

    //on résout le schéma en temps
    while t < T {
        //on détermine lambda
        let mut lambda1=vec![0.;N];
        let mut lambda2=vec![0.;N];
        for i in 1..N+1 {
            lambda1[i-1] = lambda(uapp[i],uapp[i+1], v0, pho0);
            lambda2[i-1] = lambda(uapp[i-1],uapp[i], v0, pho0);
        }
        let mut lambda = 0.;
        for i in 0..N {
            if lambda1[i]>lambda {
                lambda=lambda1[i];
            }
            if lambda2[i]>lambda {
                lambda=lambda2[i];
            }
        }

        let mut dt = dx/lambda;

        //on détermine u en n+1
        for i in 1..N+1 {
            let flux1 = flux(uapp[i],uapp[i+1], v0, pho0);
            let flux2 = flux(uapp[i-1],uapp[i], v0, pho0);
            uapp_p1[i]=uapp[i]-dt/dx*(flux1-flux2);
        }
        
        //on copie le vecteur u en n+1 déterminé avant dans uapp 
        uapp[0]=0.; //condition au bord u(a,t)=1
        for i in 1..N+1 {
            uapp[i]=uapp_p1[i];
        }
        uapp[N+1]=uapp[N];

        //est-ce que le trafic est devenu complètement fluide pho=0 ?
        //si pour tout i pho=0 alors oui :
        let mut fluide=true;
        for i in 0..N+2 {
            if uapp[i]>10e-8 {
                fluide=false;
                break;
            }
        }
        //si le trafic est fluide on peut sortir du while
        if fluide==true {
            println!("Au bout de {:?} minutes le trafic est devenu complètement fluide.",t*60.);
            break;
        }


        t=t+dt;
    }

    (x,uapp,t)
}

fn plotu(N: usize, L: f64, T: f64, v0: u32, pho0: u32) -> std::io::Result<()> {
    let fic = File::create("resu_rusanov_trafic_t2.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...

    let (xn,uapp,tn)=solve_rusanov_trafic(N,L,T,v0,pho0);

    writeln!(fic,"{:?}",tn)?;

    for i in 0..N+2 {
        writeln!(fic, "{} {}", xn[i], uapp[i])?;
    }

    Ok(())
}

fn main() {
    let L=1.;
    let v0=130;
    let pho0=200;
    let N=400;
    let T=3.; //en minutes
    plotu(N, L, T/60., v0, pho0); //on convertit T en heures
}
