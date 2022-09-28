use std::fs::File;
use std::io::BufWriter;
use std::io::{Write};

const C: f64 = 1.;

//QUESTION 4

//Solution exacte en fonction de x et de t
fn u_exact(x: f64, t: f64) -> f64 {
    if x > C*t {
        return 0.;
    } 
    else {
        return f64::exp(x/C-t);
    }
}

//Conditions aux limites 
//si lim=0, on retourne u(0,t)
//si lim=1, on retourne u(x,0)
fn conditions_limites(t: f64, lim: i32) -> f64 {
    if lim==0 { //u(0,t)
        return f64::exp(-t);
    }
    else { //u(x,0)
        return 0.;
    }
}

// [u_0^n,0,...,0]
fn G(t: f64, N: usize) -> Vec<f64> {
    let mut Gn = vec![0.;N];
    Gn[0]=conditions_limites(t,0);
    Gn
}

//produit matrice vecteur
fn produit_matvec(l: &Vec<f64>, d: &Vec<f64>, u: &Vec<f64>, x: &Vec<f64>) -> Vec<f64> {
    let n = d.len();
    let mut b = vec![0.;n];

    b[0]=d[0]*x[0]+u[0]*x[1];
    for i in 1..n-1 {
        b[i] = l[i-1]*x[i-1]+d[i]*x[i]+u[i]*x[i+1];
    }
    b[n-1] = l[n-2]*x[n-2]+d[n-1]*x[n-1];

    b
}

//Résolution de l'équation de transport
fn solve(L: f64, N: usize, T: f64, dt: f64) -> (Vec<f64>,f64) {
    let dx=L/(N as f64);
    let mut x = vec![0.;N];
    for i in 0..N {
        x[i]=dx*((i+1) as f64);
    } 
    // U0
    let mut u_app=vec![0.;N];
    let mut tn = 0.;
    // A
    let low_A = vec![-1.;N-1];
    let diag_A = vec![1.;N];
    let sup_A = vec![0.;N-1];
    // produit
    let mut b = vec![0.;N];
    let mut Gn = vec![0.;N];
    while tn < T {
        // produit AU^n
        b=produit_matvec(&low_A, &diag_A, &sup_A, &u_app);
        //Gn
        Gn=G(tn,N);
        for i in 0..N {
            u_app[i]=u_app[i] + C*dt/dx*(Gn[i]-b[i]);
        }

        tn=tn+dt;
    }
    (u_app,tn)
}

// affichage solution
fn plotu(L: f64, N: usize, T: f64, dt: f64) -> std::io::Result<()> {
    let dx=L/(N as f64);

    let (uapp,tn) = solve(L,N,T,dt);

    let fic = File::create("resu.dat")?;
    let mut fic = BufWriter::new(fic);    

    writeln!(fic,"{:?}",tn)?;

    for i in 0..N {
        let xi = dx*((i+1) as f64);

        write!(fic,"{:?} {:?} {:?}", xi , uapp[i] , u_exact(xi, T) )?;
        writeln!(fic)?;
    }

    Ok(())
}

//QUESTION 5

//calcule la norme L^1
fn norme_L1(x: Vec<f64>) -> f64{
    let n=x.len();
    let mut norme = 0.;
    for i in 0..n {
        norme = norme + f64::abs(x[i]);
    }
    norme
}

//calcule la norme L^2
fn norme_L2(x: Vec<f64>) -> f64{
    let n=x.len();
    let mut norme = 0.;
    for i in 0..n {
        norme = norme + f64::powi(f64::abs(x[i]),2);
    }
    f64::sqrt(norme)
}

//soustraction de deux vecteur
fn difference(v1: &Vec<f64>, v2: &Vec<f64>) -> Vec<f64> {
    let n=v1.len();
    let mut diff = vec![0.;n];
    for i in 0..n {
        diff[i]=v1[i]-v2[i];
    }
    diff
}

fn convergence(L: f64, tab_N: Vec<usize>, T: f64) -> std::io::Result<()> {
    let nb_N = tab_N.len();

    let fic = File::create("cvg.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...

    //On calcule les erreurs
    let mut tab_erreur_L1 = vec![0.;nb_N];
    let mut tab_erreur_L2 = vec![0.;nb_N];

    for n in 0..nb_N {
        let N=tab_N[n];
        let dx = L/(N as f64);
        let dt=0.9*dx/C;

        let (uapp,tn) = solve(L, N, T, dt);
        let mut uex = vec![0.;N];
        for i in 0..N {
            let xi = dx*((i+1) as f64);
            uex[i] = u_exact(xi, T);
        }

        tab_erreur_L1[n] = norme_L1(difference(&uapp,&uex)) * dx ;
        tab_erreur_L2[n] = norme_L2(difference(&uapp,&uex)) * dx;           
    }

    //On écrit dans le fichier
    for n in 0..nb_N {
        writeln!(fic,"{:?} {:?} {:?}",tab_N[n],tab_erreur_L1[n],tab_erreur_L2[n])?;
    }

    Ok(())
}

fn main() {
    let L=1.;
    let T = 0.5;
    let tab_N = vec![10,50,100,200,300,500,1000,2000,3000];
    convergence(L, tab_N, T);
}

// QUESTION 6

//Résolution de l'équation de transport (schéma centré)
fn schema_centre(L: f64, N: usize, T: f64, dt: f64) -> (Vec<f64>,f64) {
    let dx=L/(N as f64); //
    let mut x = vec![0.;N];
    for i in 0..N {
        x[i]=dx*((i+1) as f64);
    } 
    // U0
    let mut u_app=vec![0.;N];
    let mut tn = 0.;
    // A
    let low_A = vec![-1.;N-1];
    let diag_A = vec![0.;N];
    let sup_A = vec![1.;N-1];
    // produit
    let mut b = vec![0.;N];
    let mut Gn = vec![0.;N];
    while tn < T {
        // produit AU^n
        b=produit_matvec(&low_A, &diag_A, &sup_A, &u_app);
        //Gn
        Gn=G(tn,N);
        for i in 0..N {
            u_app[i]=u_app[i] + C*dt/2./dx*(Gn[i]-b[i]);
        }

        tn=tn+dt;
    }
    (u_app,tn)
}

fn plotu_centre(L: f64, N: usize, T: f64, dt: f64) -> std::io::Result<()> {
    let dx=L/(N as f64);

    let (uapp,tn) = schema_centre(L,N,T,dt);

    let fic = File::create("resu_centre.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...    

    writeln!(fic,"{:?}",tn)?;

    for i in 0..N {
        let xi = dx*((i+1) as f64);

        write!(fic,"{:?} {:?} {:?}", xi , uapp[i] , u_exact(xi, tn) )?;
        writeln!(fic)?;
    }

    Ok(())
}


// QUESTION 7

fn schema_lax_wendrof(L: f64, N: usize, T: f64, dt: f64) -> (Vec<f64>,f64) {
    let dx=L/(N as f64); //
    let mut x = vec![0.;N];
    for i in 0..N {
        x[i]=dx*((i+1) as f64);
    } 
    let alpha = C*dt/dx;    
    // U0
    let mut u_app=vec![0.;N];
    let mut tn = 0.;
    // A
    let low_A = vec![(alpha+1.)/2.;N-1];
    let mut diag_A = vec![-alpha;N]; diag_A[N-1]=-(alpha+1.)/2.;
    let sup_A = vec![(alpha-1.)/2.;N-1];
    // produit
    let mut b = vec![0.;N];
    let mut Gn = vec![0.;N];
    while tn < T {
        // produit AU^n
        b=produit_matvec(&low_A, &diag_A, &sup_A, &u_app);
        //Gn
        Gn=G(tn,N);

        for i in 0..N {
            u_app[i] = u_app[i] + alpha*b[i]+alpha*(alpha+1.)/2.*Gn[i];
        }

        tn=tn+dt;
    }
    (u_app,tn)
}

fn plotu_lax_wendrof(L: f64, N: usize, T: f64, dt: f64) -> std::io::Result<()> {
    let dx=L/(N as f64);

    let (uapp,tn) = schema_lax_wendrof(L,N,T,dt);

    let fic = File::create("resu_lax.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...    

    writeln!(fic,"{:?}",tn)?;

    for i in 0..N {
        let xi = dx*((i+1) as f64);

        write!(fic,"{:?} {:?} {:?}", xi , uapp[i] , u_exact(xi, tn) )?;
        writeln!(fic)?;
    }

    Ok(())
}

//Dans L^infini
fn norme_Linf(x: Vec<f64>) -> f64{
    let n=x.len();
    let mut norme = 0.; //=max(|x_i|)
    for i in 0..n {
        if f64::abs(x[i])>norme {
            norme=f64::abs(x[i])
        }
    }
    norme
}

fn plot_erreur_inf(L: f64, tab_N: Vec<usize>, T: f64) -> std::io::Result<()> {
    let nb_N = tab_N.len();

    let fic = File::create("lax_Linf.dat")?;
    let mut fic = BufWriter::new(fic); // create a buffer for faster writes...

    //On calcule l'erreur en norme infini
    let mut tab_erreur_Linf = vec![0.;nb_N];

    for n in 0..nb_N {
        let N=tab_N[n];
        let dx = L/(N as f64);
        let dt=dx/C;

        let (uapp,tn) = solve(L, N, T, dt);
        let mut uex = vec![0.;N];
        for i in 0..N {
            let xi = dx*((i+1) as f64);
            uex[i] = u_exact(xi, tn);
        }

        tab_erreur_Linf[n] = norme_Linf(difference(&uapp,&uex));       
    }

    //On écrit dans le fichier
    for n in 0..nb_N {
        writeln!(fic,"{:?} {:?}",tab_N[n],f64::ln(tab_erreur_Linf[n]))?;
    }

    Ok(())
}

// fn main() {
//     let L = 1.;
//     let T = 0.5;
//     let N = 200;
//     let dx = L/(N as f64);
//     let dt=0.5*dx/f64::abs(C); // < dx
//     plotu_lax_wendrof(L, N, T, dt);
//     let tab_N = vec![50,100,200,300,500,1000,2000,3000];
//     plot_erreur_inf(L, tab_N, T);
// }
