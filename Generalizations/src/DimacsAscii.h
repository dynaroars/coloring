static char preamble[10000];//stores DIMACS graph file preamble

void get_params(int &nVertices,int &nEdges){

		/* getting Nr_vert and Nr_edge from the preamble string, 
							 containing Dimacs format "p ??? num num" */

		char c, *tmp;
		char * pp = preamble;
		int stop = 0;
		tmp = (char *)calloc(100, sizeof(char));
		
		
		while (!stop && (c = *pp++) != '\0'){
			switch (c)
			  {
				case 'c':
				  while ((c = *pp++) != '\n' && c != '\0');
				  break;
				  
				case 'p':
				  sscanf(pp, "%s %d %d\n", tmp, &nVertices, &nEdges);
				  stop = 1;
				  break;
				  
				default:
				  break;
			  }
		}
		
		free(tmp);
}


void read_graph_DIMACS_ascii(const char *file,int &nVertices,int &nEdges, 
							 vector<vector<bool> >&edgeMatrix, 
							 vector<vector<int> >&weightMatrix,
							 vector<int> &ku){

	int c, oc;
	char * pp = preamble;
	int i=0,j=0,w=1;
	int uu,kk=1;
	FILE *fp;

	
	if ( (fp=fopen(file,"r"))==NULL )
	  { printf("ERROR: Cannot open infile\n"); exit(10); }

	for(oc = '\0' ;(c = fgetc(fp)) != EOF && (oc != '\n' || c != 'e')
		; oc = *pp++ = c);
 
	ungetc(c, fp); 
	*pp = '\0';
	get_params(nVertices,nEdges);

  //initalize the edge matrix, big waste of space
	for(int i =0;i<nVertices;++i){
	  vector<bool>temp;
	  edgeMatrix.push_back(temp);

	  vector<int>temp1;
	  weightMatrix.push_back(temp1);

	  for(int j=0;j<nVertices;++j){
		edgeMatrix.at(i).push_back(false);
		weightMatrix.at(i).push_back(-10);//default is no different, (i.e., no edge)
	  }

	  ku.push_back(-1);//default size is 1
	}

	while ((c = fgetc(fp)) != EOF){
		switch (c)
		  {
			case 'e':
			  if (!fscanf(fp, "%d %d %d", &i, &j, &w))
				{ printf("ERROR: corrupted inputfile\n"); exit(10); }
			  
				
			  edgeMatrix.at(i-1).at(j-1)=true;
			  edgeMatrix.at(j-1).at(i-1)=true;
			  
				//printf("e %d %d %d\n",i,j,w);
			  
			  weightMatrix.at(i-1).at(j-1)=w;
			  weightMatrix.at(j-1).at(i-1)=w;
			  
			  break;
			  
		  case 'n':
			if (!fscanf(fp, "%d %d", &uu, &kk)) { printf("ERROR: corrupted inputfile\n"); exit(10); }
			ku.at(uu-1)=kk;
			//printf("n %d %d\n",uu,kk);
			break;

		  case '\n':
			
		  default:
			break;
		  }
	}
	
	
/* 	for(int i = 0 ;  i< weightMatrix.size();++i){ */
/* 	  for(int j = 0 ;  j< weightMatrix.at(i).size();++j){ */
/* 		printf("w %d",weightMatrix.at(i).at(j)); */
/* 	  } */
/* 	} */

	fclose(fp);

}



