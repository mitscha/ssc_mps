# Makefile to generate fig2.pdf, fig3.pdf, fig4.pdf

all: fig2.pdf fig3.pdf fig4.pdf

fig2.pdf: CE-tvsrho_omp.dat CE-tvsrho_mp.dat CE-tvssig_omp.dat CE-tvssig_mp.dat CE-rhovssig_omp.dat CE-rhovssig_mp.dat
	pdflatex -halt-on-error -interaction=batchmode fig2.tex

fig3.pdf: TDOMP-0.20.dat TDMP-0.20.dat TDOMP-0.50.dat TDMP-0.50.dat TPROMP-0.20.dat TPRMP-0.20.dat TPROMP-0.50.dat TPRMP-0.50.dat FPROMP-0.20.dat FPRMP-0.20.dat FPROMP-0.50.dat FPRMP-0.50.dat
	pdflatex -halt-on-error -interaction=batchmode fig3.tex

fig4.pdf: CEs-itersens_faces.dat CEs-itersens_synth_0.50.dat TPFPs-itersens_faces.dat TPFPs-itersens_synth_0.50.dat TF1norm-itersens_faces.dat TF1norm-itersens_synth_0.50.dat
	pdflatex -halt-on-error -interaction=batchmode fig4.tex



CE-tvsrho_omp.dat:
	matlab -nojvm -r "phasediag(); exit(0)"
CE-tvsrho_mp.dat:
	
CE-tvssig_omp.dat:
	
CE-tvssig_mp.dat:
	
CE-rhovssig_omp.dat:
	
CE-rhovssig_mp.dat:
	


TDOMP-0.20.dat:
	matlab -nojvm -r "ROCtau; exit(0)"
TDMP-0.20.dat:
	
TDOMP-0.50.dat:
	
TDMP-0.50.dat:
	
TPROMP-0.20.dat:
	
TPRMP-0.20.dat:
	
TPROMP-0.50.dat:
	
TPRMP-0.50.dat:
	
FPROMP-0.20.dat:
	
FPRMP-0.20.dat:
	
FPROMP-0.50.dat:
	
FPRMP-0.50.dat:
	

CEs-itersens_faces.dat:
	 matlab -nojvm -r "itersensitivity(1,0); exit(0)"
CEs-itersens_synth_0.50.dat:
	matlab -nojvm -r "itersensitivity(2,0.5); exit(0)"
TPFPs-itersens_faces.dat:
	
TPFPs-itersens_synth_0.50.dat:
	
TF1norm-itersens_faces.dat:
	
TF1norm-itersens_synth_0.50.dat:
	






clean:
	rm *.log *.gz *.aux
cleanall:
	rm *.dat *.log *.pdf *.gz *.aux
