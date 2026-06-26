#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// TODO: change alilength calculation to let fasta alignments that are wrapped to be used as input (line 65)
char **allAlignedSeqs;
int alilength = 0;
int threshold;          //less OR equal than threshold variants allowed
int cutoff;             //cutoff is the smallest detectable HOR (smaller are removed, equal are kept)
int split;

bool compareAB(int *indA, int *indB, int *snvCounto);



int main(int argc, char *argv[])
{
    printf( "***********************************************************************\n"
            "*                     HOR identification software                     *\n"
            "*                                V3.0                                 *\n"
            "*        Single file input split into repeats from two regions        *\n"
            "*  Repeats need to have their direction set in the name as D1 or D2!  *\n"
            "***********************************************************************\n\n\n");
    clock_t t;
    t = clock();

    if(argc == 6)
    {
        printf("Comparing repeats from the folder: %s\nfile: %s, repeats split after repeat no: %s, threshold: %s and cutoff value: %s \n\n", argv[1], argv[2], argv[3], argv[4], argv[5]);
    }
    else
    {
        printf("Wrong amount of arguments supplied.\nProvide (in order):\nDirectory of .fasta alignments to be compared and output placed\nName of the alignment\nWhere the alignment should be split\nThreshold value\nCutoff value\n\n");
        return 1;
    }

    split = atoi(argv[3]) - 1;
    threshold = atoi(argv[4]);
    cutoff = atoi(argv[5]);

    char alifile[1000];

    strcpy(alifile, argv[1]);
    strcat(alifile, argv[2]);


    printf("Alifile: %s\n", alifile);
    // Count how many sequences and what is the length of the alignment
    printf("Checking the alignment file and calculating the number of sequences ...\n\n");

    FILE *alignment;
	alignment = fopen(alifile, "r");

	if(alignment == NULL)
    {
        printf("Error: the file does not exist, is corrupted, or denied access\n\n");
        exit(0);
    }

    for(int c = fgetc(alignment); c != 10; c = fgetc(alignment)) {   }
    for(int c = fgetc(alignment); c != '>'; c = fgetc(alignment))
    {
        if(c != 10)
        {
            alilength = alilength + 1;
        }
    }
    printf("The length of the alignment is: %i\n\n", alilength);

	int alicount = 2;
	for (int c = fgetc(alignment); c != EOF; c = fgetc(alignment))
	{
		if(c == '>')
		{
			alicount = alicount + 1;
		}
	}
	printf("Count of aligned sequences is: %i\n\n", alicount);

	fclose(alignment);



    alignment = fopen(alifile, "r");                                                                                // read aligned sequences into a new allAlignedSeqs string array

    allAlignedSeqs = malloc(alicount * sizeof(char*));
    printf("Allocating memory for the alignment ...\n\n");
    char dir[alicount];

    char temp[alilength+1];
    for(int i = 0; i < alicount; i++)
    {
        allAlignedSeqs[i] = malloc((alilength + 1) * sizeof(char));
    }
    printf("Reading the alignment ...\n\n");
	for(int j = 0; j < alicount; j ++)
	{
		int readPos = 0;
        for(char c = fgetc(alignment); c != 'D'; c = fgetc(alignment))
		{   }
        dir[j] = fgetc(alignment);
		for(char c = fgetc(alignment); c != 10; c = fgetc(alignment))
		{	}
		for(char c = fgetc(alignment); ((c != '>') && (c != EOF)); c = fgetc(alignment))
		{
			temp[readPos] = c;
			readPos++;
		}
        strcpy(allAlignedSeqs[j], temp);

	}
    fclose(alignment);

    //prepare the output file
    FILE *csvfile;                                                                                                          // Save the output
    char csvfilename[1000];
    sprintf(csvfilename,"%sHORs_%s_t_%i_c_%i.csv", argv[1], argv[2], threshold, cutoff);
	csvfile = fopen(csvfilename, "w");
    fprintf(csvfile, "start_A,end_A,start_B,end_B,direction(1=para_2=perp),total_variant\n");
	//later do:
	//fprintf(csvfile, "%i,%i,%i,%i,%i,%i,%i\n", i+1, horSTARTAdata, horENDAdata, horSTARTBdata, horENDBdata, horDIRdata,horSNVdata);

	//prepare variables to store HOR stats currently processed

    int i,j,snvcount;
    int openHOR = 0;
    int snvBack = 0;

    printf("\n\nquadrant I: ");

    snvcount = 0;
    openHOR = 0;
    for(int i = 0; i < (alicount - split)/100; i++)
    {
        printf("o");
    }
    printf("o\nquadrant I: ");

    for(int s = split+1; s < alicount ; s++)
    {
        snvcount = 0;
        i = 0;
        j = s;
        if(j%100 == 0)
        {
            printf("I");
        }
        while(j > split && i <= split)
        {
            //printf("i %i, j %i\n", i, j);
            snvBack = snvcount;
            bool isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] != dir[j]);

            if(isSimilar && openHOR == 0)
            {
                //printf("open loop");
                openHOR = 1;
            } else if(isSimilar && openHOR > 0)
            {
                openHOR = openHOR + 1;
            } else if(!isSimilar && openHOR > 0)
            {
                if(openHOR >= cutoff)
                {
                    fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR, i - 1, j + 1, j + openHOR, 2, snvBack);
                }
                openHOR = 0;
                snvcount = 0;
            }
            i++;
            j--;
        }
        if(isSimilar && openHOR >= cutoff) // close an open HOR if it was opened but not closed due to running into the edge of the array
        {
            fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR, i - 1, j + 1, j + openHOR, 2, snvBack);
        }
        openHOR = 0;
        snvcount = 0;
    }
    snvcount = 0;
    openHOR = 0;
    printf("\nquadrant II: ");
    for(int i = 0; i < split/100; i++)
    {
        printf("o");
    }
    printf("\nquadrant II: ");
    for(int s = 1; s <= split; s++)
    {
        snvcount = 0;
        i = s;
        j = alicount - 1;
        if(i%100 == 0)
        {
            printf("I");
        }

        while(j > split && i <= split)
        {
            //printf("i %i, j %i\n", i, j);
            snvBack = snvcount;
            bool isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] != dir[j]);

            if(isSimilar && openHOR == 0)
            {
                //printf("open loop");
                openHOR = 1;
            } else if(isSimilar && openHOR > 0)
            {
                openHOR = openHOR + 1;
            } else if(!isSimilar && openHOR > 0)
            {
                if(openHOR >= cutoff)
                {
                    fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR, i - 1, j + 1, j + openHOR, 2, snvBack);
                }
                openHOR = 0;
                snvcount = 0;
            }
            i++;
            j--;
        }
        if(isSimilar && openHOR >= cutoff) // close an open HOR if it was opened but not closed due to running into the edge of the array
        {
            fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR, i - 1, j + 1, j + openHOR, 2, snvBack);
        }
        openHOR = 0;
        snvcount = 0;
    }

    printf("\nquadrant III: ");
    for(int i = 0; i < (alicount - split)/100; i++)
    {
        printf("o");
    }
    printf("o\nquadrant III: ");
    snvcount = 0;
    openHOR = 0;
    for(int s = alicount - 1; s > split; s--)
    {
        snvcount = 0;
        i = 0;
        j = s;
        if(j%100 == 0)
        {
            printf("I");
        }
        while(j < alicount && i <= split)
        {
            //printf("i %i, j %i\n", i, j);
            snvBack = snvcount;
            bool isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] == dir[j]);

            if(isSimilar && openHOR == 0)
            {
                //printf("open loop");
                openHOR = 1;
            } else if(isSimilar && openHOR > 0)
            {
                openHOR = openHOR + 1;
            } else if(!isSimilar && openHOR > 0)
            {
                if(openHOR >= cutoff)
                {
                    fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR, i - 1, j - openHOR, j - 1, 1, snvBack);
                }
                openHOR = 0;
                snvcount = 0;
            }
            i++;
            j++;
        }
        if(isSimilar && openHOR >= cutoff)
        {
            fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR, i - 1, j - openHOR, j - 1, 1, snvBack);
        }
        openHOR = 0;
        snvcount = 0;
    }
    snvcount = 0;
    openHOR = 0;
    printf("\nquadrant IV: ");
    for(int i = 0; i < split/100; i++)
    {
        printf("o");
    }
    printf("\nquadrant IV: ");
    for(int s = 1; s <= split; s++)
    {
        snvcount = 0;
        i = s;
        j = split + 1;
         if(i%100 == 0)
        {
            printf("I");
        }
        while(j < alicount && i <= split)
        {
            //printf("i %i, j %i\n", i, j);
            snvBack = snvcount;
            bool isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] == dir[j]);

            if(isSimilar && openHOR == 0)
            {
                //printf("open loop");
                openHOR = 1;
            } else if(isSimilar && openHOR > 0)
            {
                openHOR = openHOR + 1;
            } else if(!isSimilar && openHOR > 0)
            {
                if(openHOR >= cutoff)
                {
                    fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR, i - 1, j - openHOR, j - 1, 1, snvBack);
                }
                openHOR = 0;
                snvcount = 0;
            }
            i++;
            j++;
        }
        if(isSimilar && openHOR >= cutoff)
        {
            fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR, i - 1, j - openHOR, j - 1, 1, snvBack);
        }
        openHOR = 0;
        snvcount = 0;
    }
    printf("\n");

   for(int i = 0; i < alicount; i++)
    {
        free(allAlignedSeqs[i]);
    }
    free(allAlignedSeqs);

    printf("\n\nHOR output file created:\n\n%sHORs_%s_t_%i_c_%i.csv\n\n", argv[1], argv[2], threshold, cutoff);
    fclose(csvfile);
    t = clock() - t;

    FILE *log;                                                                                                          // Save the output
    char logfile[1000];
    sprintf(logfile,"%slog_%s_t_%i_c_%i.txt", argv[1], argv[2], threshold, cutoff);
	log = fopen(logfile, "w");
    fprintf(log, "HOR output file created:\n\n%sHORs_%s_t_%i_c_%i.csv\n\n", argv[1], argv[2], threshold, cutoff);
    fprintf(log, "The program finished in %f minutes.\n",((float)t)/CLOCKS_PER_SEC/60);
    fclose(log);
    printf ("The program finished in %f minutes.\n",((float)t)/CLOCKS_PER_SEC/60);

	return 0;

}





bool compareAB(int *indA, int *indB, int *snvCounto)
{
    int snv = 0;
    for(int i = 0; i < alilength; i++)
    {
        snv = snv + (allAlignedSeqs[*indA][i] != allAlignedSeqs[*indB][i]);
    }
    if(snv <= threshold)
    {
        *snvCounto = *snvCounto + snv;
    }
    //printf ("\n comparing %i with %i,  snc count is %i", *indA, *indB, snv);
    return(snv <= threshold);
}
