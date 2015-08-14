/* bb.c */
double planck(double t, double freqmin, double freqmax);
double get_rand_pow(double x1, double x2, double alpha);
double get_rand_exp(double alpha_min, double alpha_max);
double integ_planck_d(double alphamin, double alphamax);
int init_integ_planck_d(void);
double planck_d(double alpha);
double emittance_bb(double freqmin, double freqmax, double t);
double check_fmax(double fmin, double fmax, double temp);
/* get_atomicdata.c */
int get_atomic_data(char masterfile[]);
int index_lines(void);
int index_phot_top(void);
int index_inner_cross(void);
int index_collisions(void);
void indexx(int n, float arrin[], int indx[]);
int limit_lines(double freqmin, double freqmax);
int check_xsections(void);
/* python.c */
int main(int argc, char *argv[]);
int help(void);
int init_geo(void);
int photon_checks(PhotPtr p, double freqmin, double freqmax, char *comment);
int get_spectype(int yesno, char *question, int *spectype);
int qdisk_init(void);
int qdisk_save(char *diskfile, double ztot);
int read_non_standard_disk_profile(char *tprofile);
int init_advanced_modes(void);
/* photon2d.c */
int translate(WindPtr w, PhotPtr pp, double tau_scat, double *tau, int *nres);
int translate_in_space(PhotPtr pp);
double ds_to_wind(PhotPtr pp);
int translate_in_wind(WindPtr w, PhotPtr p, double tau_scat, double *tau, int *nres);
int walls(PhotPtr p, PhotPtr pold);
/* photon_gen.c */
int define_phot(PhotPtr p, double f1, double f2, long nphot_tot, int ioniz_or_final, int iwind, int freq_sampling);
double populate_bands(double f1, double f2, int ioniz_or_final, int iwind, struct xbands *band);
int xdefine_phot(double f1, double f2, int ioniz_or_final, int iwind);
int xmake_phot(PhotPtr p, double f1, double f2, int ioniz_or_final, int iwind, double weight, int iphot_start, int nphotons);
double star_init(double r, double tstar, double freqmin, double freqmax, int ioniz_or_final, double *f);
int photo_gen_star(PhotPtr p, double r, double t, double weight, double f1, double f2, int spectype, int istart, int nphot);
double disk_init(double rmin, double rmax, double m, double mdot, double freqmin, double freqmax, int ioniz_or_final, double *ftot);
int photo_gen_disk(PhotPtr p, double weight, double f1, double f2, int spectype, int istart, int nphot);
int phot_gen_sum(char filename[], char mode[]);
double bl_init(double lum_bl, double t_bl, double freqmin, double freqmax, int ioniz_or_final, double *f);
/* saha.c */
int nebular_concentrations(PlasmaPtr xplasma, int mode);
int concentrations(PlasmaPtr xplasma, int mode);
int saha(PlasmaPtr xplasma, double ne, double t);
int lucy(PlasmaPtr xplasma);
int lucy_mazzali1(double nh, double t_r, double t_e, double www, int nelem, double ne, double density[], double xne, double newden[]);
int fix_concentrations(PlasmaPtr xplasma, int mode);
double get_ne(double density[]);
/* spectra.c */
int spectrum_init(double f1, double f2, int nangle, double angle[], double phase[], int scat_select[], int top_bot_select[], int select_extract, double rho_select[], double z_select[], double az_select[], double r_select[]);
int spectrum_create(PhotPtr p, double f1, double f2, int nangle, int select_extract);
int spectrum_summary(char filename[], char mode[], int nspecmin, int nspecmax, int select_spectype, double renorm, int loglin);
int spectrum_restart_renormalise(int nangle);
/* wind2d.c */
int define_wind(void);
int where_in_grid(double x[]);
int vwind_xyz(PhotPtr p, double v[]);
int wind_div_v(WindPtr w);
double rho(WindPtr w, double x[]);
int mdot_wind(WindPtr w, double z, double rmax);
int get_random_location(int n, int icomp, double x[]);
int zero_scatters(void);
int check_corners_inwind(int n);
int set_nstart_nstop(void);
/* wind.c */
int where_in_wind(double x[]);
int wind_check(WindPtr www, int n);
double model_velocity(int ndom, double x[], double v[]);
int model_vgrad(int ndom, double x[], double v_grad[][3]);
double model_rho(int ndom, double x[]);
/* vector.c */
double dot(double a[], double b[]);
double length(double a[]);
int renorm(double a[], double scalar);
int cross(double a[], double b[], double 