<?php
/*
This code will be used to insert E. coli strains from Nichols et al, 2011 into the website
Execute by: php NicholsStrainCreations.php -w /Library/WebServer/Documents/omp/peter/

The test website's URL: https://microbialphenotypes.org/peter/index.php/Main_Page
*/
$params = getopt( "w:" );
set_IP($params['w']);
$maintClass = "NicholsMaint";
class NicholsMaint extends Maintenance {
    public function __construct() {
        parent::__construct();
        $this->parse_parameters();
    }

    public function execute() {
        $filename = 'JW123.txt';
        echo "opening $filename\n";
        $fh = fopen($filename, 'r');


        if(!$fh) die;


        $pageNo=0; #This is a counter. No need to change it
        $totalPageNo=3955; //change this to limit the No. of pages being created  #Change to 3818 later
        while ($line = fgets($fh)) {


            #(Vefify with Sandy) If there is an empty line, ignore it (but why do I need trim($line) == ""?)
            if ($line{0} == "/n" || trim($line) == "") {
                continue;
            }


            $trimmedLine=trim($line);


            #this makes variables for all of the data from the file
            list ($strain, $ECK, $gene, $allele, $bkd, $category, $pubmedID, $JWSynonym, $CGSCurl) = explode("\t", $trimmedLine);




            #This correct $strain ("BW25113 + allele in column D" is replaced with real allele)
            if(preg_match('/\+ allele in column D/',$strain)){ #Should use '\+' instead of just '+'. I don't know why preg_match needs '/' at the begining and the end of the string
               $strain=str_replace('+ allele in column D',$allele, $strain);
            }



            #Have to re-define synonyms
            if($JWSynonym!="NA"){
                $synonyms=$JWSynonym;
            }else{
                $synonyms="";
            }




#==============================This part checks if the strain page(s) already exists. If so, skip=======================

            #I felt reluctant to use this function (deprecated in PHP5.5.0) but couldn't find an alternative
            ## Note that there are ' in some allele names (which should be escaped). But they all have JW numbers (' won't appear in page names) => doesn't matter here so I will comment out that function
            ##
            #$sql_strain=mysql_real_escape_string($strain);
            ## This is to test if any special characters in the alleles would cause problems if they are not escaped: echo mysql_real_escape_string("abc def ' - ( )  : :  + -")."\n";


            $sql_strain=str_replace(" ","_",$strain); #Even if escaped by mysql_real_escape_string(), the " " still won't change to "_" => Have to do it manually


            $dbr = wfGetDB( DB_SLAVE );
            $result = $dbr->select(
                array('page'),
                '*',
                array("page_title LIKE '%Escherichia_coli_K-12_".$sql_strain."%'"),
                __METHOD__,
                array()
            );


            $Exist=0;
            #Note: If JWXXXX-Y has been redirected to just JWXXXX, there will be 2 pages found. However, the JWXXXX-Y page will have no content
            foreach($result as $x) { #If no pages are found, this foreach() is not run. I don't know why


                #Uses the pages found and creates a title object and finds the tables on the strain pages that were previously created
                $newtitle = Title::newFromText($x->page_title);
                $new_page = new WikiPageTE($newtitle);

                echo $newtitle."\n";
                echo "->This page exists and will not be made"."\n";
                $Exist=1; #If the page is found, set this variable to 1

            }

            if($Exist==1) continue; #If the page(s) already exists, skip to the next while() iteration

#=======================================================================================================================


            #creates a new id number for each strain
            $dbw = wfGetDB(DB_SLAVE);
            $dbw->insert('omp_master.strain', array('id' => null, 'name' => $strain));  #omp_master is a database. strain is a table in it. Here it inserts a strain named $JW  #If id is null, the number is going to be auto-incremented
            $id = $dbw->insertId(); #This pulls out the id created by the above line

            #creates a new page based on a template and names the page
            $new_page_pageName = "OMP_ST:" . "$id" . "_!_" . "Escherichia coli K-12 ".$strain; #Have to check if there are like ' in the genes_symbol that can possibly causes crashes
            $newpagetemplateTitle = Title::newFromText("Template:StrainPage");
            $templatePage = new WikiPageTE($newpagetemplateTitle);

            #gets the content for the template page
            $text1 = $templatePage->getContent();
            $reason = 'Page creation by ' . __CLASS__;

            #creates a url to link to later
            $t = Title::newFromText($new_page_pageName);
            $t->getFullURL();

            #creates the new page and saves it
            $wikiPage = new WikiPageTE($t);
            $wikiPage->save($text1, $reason,EDIT_FORCE_BOT);
            ##EDIT_FORCE_BOT is an optional argument to the save function. The purpose of having it is to not flush the recent changes
            ##Still, by going through WikiPageTE.php and a couple of related files, I wasn't able to figure out why EDIT_FORCE_BOT can make the page edits flagged by bots. I tested it and it worked, so I used it
            $wikiPage->touch(); #Undefined index:ancestry and Undefined index: genotype are because of this (Confirmed by putting die; in front of and after this line)

            #Output the created pages to pagesToBeDeleted.txt => This file can be used to deleted the crated pages as follows (By using the Mediawiki built-in deleteBatch.php):
            #php /Library/WebServer/Documents/omp/peter/maintenance/deleteBatch.php -u root -r "delete" pagesToBeDeleted.txt
            file_put_contents("pagesToBeDeleted.txt",$new_page_pageName."\n",FILE_APPEND);



            if($bkd=='parent: [[OMP ST:800 ! Escherichia coli K-12 BW25113]]') $parentname='OMP ST:800 ! Escherichia coli K-12 BW25113';
            if($bkd=='parent: [[OMP ST:159 ! Escherichia coli K-12 MG1655]]') $parentname='OMP ST:159 ! Escherichia coli K-12 MG1655';

            #uses the parent bioneer strain as the parent and gets the genotype
            $p = Title::newFromText($parentname);

            #gets strain info table from parent page for later insertion
            $temp = "Strain_info_table";
            $parent_page = new WikiPageTE($p);
            $parenttable = $parent_page->getTable($temp);
            $box = $parenttable[0];
            $parentgenotyperow = $box->get_row_hash(0);
            $wikiPage->touch();



            if($CGSCurl!="NA"){
                $strainRef='*'.$pubmedID."\n".'*'.'['.$CGSCurl.' '.'Coli Genetic Stock Center'.']';
                $strainAvail='['.$CGSCurl.' '.'Coli Genetic Stock Center'.']';
            }else {
                $strainRef=$pubmedID;
                $strainAvail=''; #(Verified) '' means blank for wiki pages
            }





            #new row on strain info table and inserts strain name, parent genotype, and parent name
            $newtable = $wikiPage->getTable("Strain_info_table");
            $box = $newtable[0];
            $box->insert_row(
                "Escherichia coli K-12 "."$strain"."||". #Strain name
                $synonyms."||".
                "*Pangenome: Escherichia coli\n*Subspecies and/or strain: ".'K-12'."\n*NCBI Taxonomy ID: ".'[http://www.ncbi.nlm.nih.gov/taxonomy?term=83333 83333]'."||". #Taxon Information
                $parentgenotyperow['genotype']." ".$allele."||". #Genotype
                $strainRef."||". #Strain Reference
                $strainAvail."||". #Strain availability
                $bkd."||" #Ancestry
            );
            $wikiPage->touch();
            echo $new_page_pageName." is made!\n";
            echo "genotype= ".$parentgenotyperow['genotype']." ".$allele."\n";


            $pageNo++;

            if ($pageNo==$totalPageNo) {
                break;
            }
        }
    }







    private function parse_parameters(){
        $this->addOption( "people", "say hello to this", $required = false, $withArg = true, $shortName = 'l' );
    }
}
require_once( RUN_MAINTENANCE_IF_MAIN );
function set_IP($path){
    global $IP;
    if ( isset($path) && is_file("$path/maintenance/Maintenance.php") ){
        $IP = $path;
        require_once( $IP . "/maintenance/Maintenance.php" );
        return $path;
    } else {
        die ("need -w <path-to-wiki-directory>");
    }
}