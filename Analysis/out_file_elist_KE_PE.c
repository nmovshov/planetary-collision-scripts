#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

/* little program that reads particle info from the SPHERAL output files that 
* Naor Movshovitz has formatted and classifies them by equivalence class
* per Numerical recipes */  

/* cf this note at http://www.geomview.org/docs/sgarchive/msg01401.shtml */

/* 12/02/13: edited out superfluous geomview stuff */

/* simple equivalence-class program based on particle distances */

/* 12/02/13: start to improve process: first largest group */

void dist(double x1, double y1, double z1, double x2, double y2, double z2,
	double *d);

main(int argc, char *argv[]) {

	int ihead, iline, nhead, icount, iflag, iter;
	int i, j, n, m, npt, file_len;
	int kk, jj, idq, inew, ngroup;
	int ngmax, nngmax, nng, ngbound, ngbound0;

	double x1, y1, z1, x2, y2, z2, distq, hfac, hh;

	double gmass, gmass0, gboundmass, xcm, ycm, zcm, vxcm, vycm, vzcm;
	double xx, yy, zz, vxx, vyy, vzz, rrmn, etot;

	double Gbig, MTarget;
	double Gphitotal, Mtotal;

	char file[80];
        char line[200];
	char string1[132], string2[132];
        FILE *fp, *out;

	char begin_head_string[3];
	char end_head_string[3];

       	if (argc <= 1) {
                puts("enter filename");
                return 1;
        }
        sscanf(argv[1], "%s",&file);

/* set some expected strings */
	strcpy(begin_head_string,"###");
	strcpy(end_head_string,"###"); 

        fp = fopen(file, "r");

/* the format seems to have 26 header lines */
	nhead = 26;

/* start to read header lines */
	ihead = 1;
	iline = 0;

	while (iline < nhead) {
		fgets(line,200,fp);
/*		printf("line read: %s",line);  */
		if (iline == 1) {
			ihead = strncmp(line,begin_head_string,3);
			if (ihead == 0) {
			printf("begin header %s detected on line : %d\n",
				line,iline); 
			}

		}
		if (iline == nhead-1 ) {
			ihead = strncmp(line,end_head_string,3);
			if (ihead == 0) {
			printf("end header %s detected on line : %d\n",
				line,iline); 
			}
		}
		iline++;

	}

	iline = 0;

/* count lines */
	while(fgets(line,200,fp) != NULL) {
		iline++;
	}
        fclose(fp);
	printf("number of data lines = %d\n",iline);

	npt = iline;

/* allocate particle data */
	int id[npt], id_eos[npt];
	double x[npt], y[npt], z[npt];
	double vx[npt], vy[npt], vz[npt];
	double mass[npt], rho[npt], p[npt], T[npt], U[npt];
	double hmin[npt], hmax[npt];

/* particle data for members in the biggest group */
	double phig[npt], zkeg[npt];

/* tag members for iterations */
	int itag[npt];

/* reopen and re-read file */
       fp = fopen(file, "r");

/* the format seems to have 26 header lines */

/* start to read header lines */
        ihead = 1;
        iline = 0;

        while (iline < nhead) {
                fgets(line,200,fp);
/*              printf("line read: %s",line);  */
                if (iline == 1) {
                        ihead = strncmp(line,begin_head_string,3);
                        if (ihead == 0) {
                        printf("begin header %s detected on line : %d\n",
                                line,iline);
                        }

                }
                if (iline == nhead-1 ) {
                        ihead = strncmp(line,end_head_string,3);
                        if (ihead == 0) {
                        printf("end header %s detected on line : %d\n",
                                line,iline);
                        }
                }
                iline++;

        }

        i = 0;
/* read lines */
        while(fgets(line,200,fp) != NULL) {
		sscanf(line,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                    %lf %lf %lf %lf %lf", &id[i],&id_eos[i],&x[i],&y[i],&z[i],
                        &vx[i],&vy[i],&vz[i],&mass[i],&rho[i],&p[i],&T[i],&U[i],
			&hmin[i],&hmax[i]);
		i++;
        }
        fclose(fp);
        printf("number of data lines = %d\n",i);

	file_len = strlen(file);
	printf("length of file name is %d\n",file_len);

/* implement equivalence list stuff from NR */
	int nf[npt], nfkey[npt], nfgroup[npt], nngroup[npt];
/* equivalence class routine from NR 3rd edition, p 440-441 */
	nf[0] = 0;
/* set dummy value? */
	for (jj=1;jj<=npt-1;jj++) {
		nf[jj] = -jj;
/*                printf("dummy nf = %d \n",nf[jj]); */
	}  /* end loop */

/* calculate target mass (id = 0) */
        MTarget = 0.0;
        for (jj=0;jj<=npt-1;jj++) {
                if (id[jj] == 0) { MTarget = MTarget + mass[jj]; }
        }  /* end loop */
        printf("Target mass = %12.4e\n",MTarget);

/* factor for distance measure maybe should be 1? less than1 ? */
	hfac = 1.0;
	printf("NB!  hfac set to %9.4f\n",hfac);

/* loop over 1st element */
	for (jj=1;jj<=npt-1;jj++) {
		if ((jj/1000)*1000 - jj == 0) { 
			printf("outer loop particle jj = %d\n",jj) ;
                        fflush(stdout);
		}
		nf[jj] = jj;
/* loop over 2nd element */
		for (kk=0;kk<jj;kk++) {
		nf[kk] = nf[nf[kk]];

/* calculate distance between points jj+1 and kk+1 */
		x1 = x[jj+1];
		y1 = y[jj+1];
		z1 = z[jj+1];

		x2 = x[kk+1];
		y2 = y[kk+1];
		z2 = z[kk+1];

		hh = 0.25 * (hmin[jj+1] + hmax[jj+1] + 
			hmin[kk+1] + hmax[kk+1]);
/* reduce the distance some */
		hh = hfac * hh;

		dist(x1,y1,z1,x2,y2,z2,&distq);

/* the idea is, if the dist between particles jj+1 and kk+1  is less than
* the avg of their smoothing lengths, put them in the same class (set 
* idq = 1) */
		idq = 0;
		if (distq < hh) idq = 1;

/* this is the magic line according to NR */
/* apply the magic line */
		if (idq == 1) nf[nf[nf[kk]]] = jj; 	

		} /*end kk loop */
	} /* end jj loop */

/* another magic line */
	for (jj=0;jj<npt;jj++) nf[jj]=nf[nf[jj]];

/*	for (i=0;i<=npt-1;i++)  {
		printf("point %d nf = %d\n",i,nf[i]);
	} */

/* At this point, the particles have allegedly been sorted. However, the
* index nf apparently lists particles according to the number of the 
* last particle in the group. Better to make a key with a group number and
* list that somehow */

/* will do this stupidly */
	ngroup = 1;
/* save particle (group) number of 1st particle */
	nfkey[ngroup] = nf[0];
	for (i=1;i<=npt-1;i++)  {
/* check to see if group # of particle has changed */
		if (nf[i] !=  nf[i-1]) {
/* set flag for possible new group */
			inew = 1;
/* if it the group number has changed, is it the same as an earlier number? */
			j = 0;
			while (j < i && inew == 1) {
				if (nf[j] == nf[i]) {
					inew = 0;
				} /* end nf comparison */
				j++;
			} /* end j while */
/* if not, create new group */
			if (inew == 1) {
				ngroup++;
				printf("new group ng=%d  i= %d nf= %d\n",
					ngroup,i,nf[i]);
				nfkey[ngroup] = nf[i];
			} /* end inew if */
		} /* end nf if */
	} /* end i loop */

/* print key number for each group */
	for (i=1;i<=ngroup;i++)  {
/*		printf(" i = %d group key %d\n",i,nfkey[i]); */
		nngroup[i] = 0;
	}

/* reassign group numbers */
	for (i=0;i<=npt-1;i++) {
		jj = nf[i];
		for (j=1;j<=ngroup;j++)  {
			if (nfkey[j] == jj) { 
				nfgroup[i] = j; 
				nngroup[j]++;
			}
		} /* end j loop */
	} /* end i loop */

/* print key number and # of particles for each group */
/*	for (i=1;i<=ngroup;i++)  {
		printf(" i = %d group key %d %d\n",i,nfkey[i],nngroup[i]);
	} */

/* number of groups */
	printf("No. of groups found = %d\n",ngroup);
/* find number of particles in each group */
	int npergroup[ngroup+1];
	npergroup[0] = 0;
	for (i=1;i<=ngroup;i++)  {
		npergroup[i] = 0;
	}
	for (n=0;n<=npt-1;n++) { 
		for (i=0;i<=ngroup;i++)  {
			if (nfgroup[n] == i) {
				npergroup[i]++;
			} /* end if */
		} /* end i loop */
	} /* end n loop */

/* print out groups with > n members */
/* and also find which group has the largest number of members */
	ngmax = 0;
	nng = 0;
	for (i=1;i<=ngroup;i++)  {
		if (npergroup[i] > nng) {
			nng = npergroup[i];
			ngmax = i;
		} /* end if */
	} /* end loop */

	for (i=1;i<=ngroup;i++)  {
		if (npergroup[i] > 10) {
			printf("group, n = %d %d\n",i,npergroup[i]);
		} /* end if */
	} /* end loop */

	printf("group %d is largest, has %d members\n",ngmax,nng);

/* ngmax is the number of the group with the most members */

/* for first iteration, tag all members of the group ( = 1), will reset to
* zero all members with positive energy for later iterations */
	for (n=0;n<=npt-1;n++) { 
		itag[n] = -1;
		if (nfgroup[n] == ngmax) {
			itag[n] = 1;
		} /* end if */
	} /* end loop */

	iflag = 1;	
	iter = 0;
	ngbound0 = nng;
	ngbound = 0;

/* put in while here to iterate */
	while (ngbound != ngbound0) {

/* find center of mass and velocity of the ngmax group */
	gmass = 0.0;
	for (n=0;n<=npt-1;n++) { 
		if (nfgroup[n] == ngmax && itag[n] == 1) {
			gmass = gmass + mass[n];
			xcm = xcm + mass[n] * x[n];
			ycm = xcm + mass[n] * x[n];
			zcm = xcm + mass[n] * x[n];

			vxcm = vxcm + mass[n] * vx[n];
			vycm = vycm + mass[n] * vy[n];
			vzcm = vzcm + mass[n] * vz[n];
		} /* end if */
	} /* end loop */

	if (iter == 0) { gmass0 = gmass; } 
	
	xcm = xcm / gmass;
	ycm = ycm / gmass;
	zcm = zcm / gmass;

	vxcm = vxcm / gmass;
	vycm = vycm / gmass;
	vzcm = vzcm / gmass;

	printf("iteration = %d =============\n",iter);
	printf("mass of group =       %12.4e\n",gmass);
	printf("c-of-m of group =     %12.4e%12.4e%12.4e\n",xcm,ycm,zcm);
	printf("vel c-of-m of group = %12.4e%12.4e%12.4e\n",vxcm,vycm,vzcm);

/* start calculating energy of group; whether members are bound */
/* MKS units */
	Gbig = 6.672e-11;

/* KE relative to group c-of-m velocity */
/* Grav PE  of group members */
	for (n=0;n<=npt-1;n++) { 

		if ((n/10000)*10000 - n == 0) { 
			printf("energy outer loop particle n = %d\n",n) ;
                        fflush(stdout);
		}

		zkeg[n] = 0.;
		phig[n] = 0.0;
		if (nfgroup[n] == ngmax && itag[n] == 1) {
			vxx = vx[n] - vxcm;
			vyy = vy[n] - vycm;
			vzz = vz[n] - vzcm;
			zkeg[n] = 0.5 * mass[n] * (vxx * vxx + vyy * vyy 
				+ vzz * vzz);
			for (m=0;m<=npt-1;m++) { 
/* Grav PE  contrib of group member m to group member n */
/* watch out not to calculate for m = n --> gives -infinity! */
		if (nfgroup[m] == ngmax &&  itag[m] == 1 && m != n) {
					xx = x[n] - x[m];
					yy = y[n] - y[m];
					zz = z[n] - z[m];
					rrmn = xx * xx + yy * yy + zz * zz;
					rrmn = sqrt(rrmn);
					phig[n] = phig[n] + mass[m] / rrmn;

				} /* end if */
			} /* end m loop */
			phig[n] = - Gbig * mass[n] * phig[n];
			
		} /* end if */
	} /* end n loop */

/* Count group members that have Grav PE + KE < 0 */	
	ngbound0 = ngbound;
       	ngbound = 0;
	gboundmass = 0.0;
	for (n=0;n<=npt-1;n++) { 
		if (nfgroup[n] == ngmax) {
			etot = phig[n] + zkeg[n];
			if (etot < 0.0 ) {
/* if particles are bound, increment counters and bound mass */
				ngbound++;
				gboundmass = gboundmass + mass[n];
			} else {
/* untag unbound particles for next iter */

/*	printf("particle %d lost; phig, ke = %12.4e%12.4e\n",
			n,phig[n],zkeg[n]); */

				itag[n] = 0;
			} /* end if */
		}  /* end if */
	} /* end loop */
	printf(" ngbound =          %d bound mass = %12.4e\n",
		ngbound,gboundmass);
	printf("number in group =   %d group mass = %12.4e\n",nng,gmass0);
	printf("fraction of bound mass = %9.4f\n",gboundmass/gmass0);
	printf("bound mass/MTarget     = %9.4f\n",gboundmass/MTarget);

/* when the number of bound particles stops changing, unset the flag */
	if (ngbound == ngbound0) {iflag = 0; }

	iter++;

	fflush(stdout);

	} /* end while */

	Gphitotal = 0.0;
	for (n=0;n<=npt-1;n++) { 
/* output position and equiv. class */
/* only for largest group that remains bound */
		if (nfgroup[n] == ngmax && itag[n] == 1) { 
			Gphitotal = Gphitotal + 0.5 * phig[n];
		} /* end if */
	} /* end loop */

	printf("Gphitotal = %12.4e\n",Gphitotal);
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("final: fraction of bound mass = %9.4f\n",gboundmass/gmass0);
        printf("final: bound mass/MTarget     = %9.4f\n",gboundmass/MTarget);

/*==================================================================== */
/* output section of program - geomview output stuff deleted */

/* simple output for clustering calculation */
	char outfile[80];
	strcpy(outfile,file);
	strcat(outfile,".out");
	printf("outfile name is %s\n",outfile);  

        out = fopen(outfile, "w"); 
	for (n=0;n<=npt-1;n++) { 
/* output position and equiv. class */
/* only for largest group that remains bound */
		if (nfgroup[n] == ngmax && itag[n] == 1) { 
			fprintf(out,"%d %13.5e%13.5e%13.5e %d %d",
				n, x[n],y[n],z[n],nf[n],nfgroup[n]);

/* add KE, PE */
			fprintf(out," %12.4e%12.4e",zkeg[n],phig[n]); 

			fprintf(out,"\n");

		} /* end if */
	} /* end loop */ 

        fclose(fp); 


}

void dist(double x1, double y1, double z1, double x2, double y2, double z2,
        double *d) {

	double dd;

	dd = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + 
		(z1 - z2) * (z1 - z2);

	*d = sqrt(dd);

}


