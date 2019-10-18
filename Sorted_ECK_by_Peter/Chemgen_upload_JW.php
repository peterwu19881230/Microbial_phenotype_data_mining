<?php

/*Goal: This file is modified from Chemgen_upload_v2.0. The difference is that this file change the ECK names to JW names. (Note: ECKs without verified JWs won't be created). 
codes for inserting strain info like this page was also written: http://microbialphenotypes.org/wiki/index.php/OMP_ST:4399_!_Schizosaccharomyces_pombe_mug20delta   */ 


//Set path to the wiki maintenance folder so we can use stuffs from there.
$params = getopt( "w:" );
set_IP($params['w']);



$maintClass = "DemoMaint";


//Class DemoMaint appears at OMPwiki by showing "Page creation by DemoMaint", so I suppose It helps create the pages.
class DemoMaint extends Maintenance {

	public function __construct() {
		parent::__construct();
		$this->parse_parameters();
	}

	
	public function execute() {
	
	
		# assign entered conditions to variables
		$strain = $this->getOption("strain","",$required = false, $withArg = true,'$strain');
		$condition = $this->getOption("condition","",$required = false, $withArg = true,'$condition');		
		$ScoreLimit = $this->getOption("ScoreLimit","",$required = false, $withArg = true,'$ScoreLimit');
		
		//Set the default value of $ScoreLimit to 0 if no value was specified
		if($ScoreLimit==NULL){
		$ScoreLimit=0;
		}
		
		#verify the content of the entered values
		echo "-------------------------------------\n";
        echo "The strain name should match: ".$strain."\n";
		echo "The condition should match: ".$condition."\n";		
		echo "The Score Limit should match: ".$ScoreLimit."\n"; 
		echo "-------------------------------------\n";
		
		
		
	
		
		
		
		// Do a database query. Conditions specified by entries.(If not specified, the result will be all. But scorelimit must be specified)
		$dbr = wfGetDB( DB_SLAVE );
		$result2 = $dbr->select(
			'chemgen.straincond',
			'*',			
			#extract data that match entered criteria
			array("cond LIKE '$condition%'","strain LIKE '%$strain%'","ABS(score)>=$ScoreLimit"),
			__METHOD__,
			array()
		);		
		
			
				
		//Ask the database to get the page title (that contains JW). 
		$dbr = wfGetDB( DB_SLAVE );
				$result1 = $dbr->select
				(					
					array('page'),
					'*',
					
					array("page_title LIKE '%JW%'"),
					__METHOD__,
					array()		
				);	
				
		
		//Parse $result1 (an object) into an array with all page titles		
		foreach($result1 as $y)
		{
		$pageTitle[]=$y->page_title;
		}			
			
	
		
		//Parse the db result into an array with only strain name strings
		foreach ($result2 as $x2)
		{
		$strain[]=$x2->strain;
		}
		//Returns a new array with only unique elements 
		$strainUnique=array_unique($strain);
		
		//change ECK names to JW if JW names exist and were verified.
		$JWarray = array_map('str_getcsv', file('Keio_strains_with_verified_JW.csv'));
		Foreach ($JWarray as $JW){
		
		//In Keio_strains_with_verified_JW.csv, $JW[0] is the ECK number, $JW[1] is the mutated gene, $JW[2] is the JW number. 
		$JW[0]=preg_quote($JW[0]);	
			
			if(preg_grep("($JW[0])",$strainUnique)!=NULL){
			
			$JWtobeCreated[]=$JW[2];
			
			
			}
		}
		print_r($JWtobeCreated);


				
		
		//$i for strains not found, $j for strains found
		$i=0; $j=0;		
		//Parse $JWtobeCreated (an object) into an array with all strain names	
		foreach ($JWtobeCreated as $x3)
		{
		
			//Because SQL will put an underscore under every space, so I have to replace strain names with space to underscore:
			$xx=str_replace(' ','_',$x3);		
			
			
			//These 2 lines escapes the modifier(they are not modifiers, but they will be interpreted as they are)?
			//$xx=str_replace('(','(\\',$xx);
			//$xx=str_replace(')',')\\',$xx);
			
			//This line quotes regular expression characters
			$xx=preg_quote($xx);			
			
			//Somehow "(" and ")" cannot work for strains that contain "()". % does the same.
			if(!preg_grep("($xx)",$pageTitle))
			{
			echo "~~~~~~~~~~~~~~~~~~~~~$xx was not found!\n";
			$i++;
			
			echo "Last page title: ".$y->page_title." -> ".$x3."  making new pages!\n";				
				
				//The following code up until "$wikiPage->touch();" is to make new pages of strains that don't exist.
				$dbw = wfGetDB( DB_SLAVE );	
				$dbw->insert( 'omp_master.strain', array( 'id' => null, 'name' => "$x3" ) );
				$id = $dbw->insertId();
	
				#creates a new page based on a template and names the page and adds the category box
				$new_page_template = "Template:StrainPage";
				$new_page_pageName = "OMP_ST:"."$id"."_!_"."Escherichia coli K-12_"."$x3";
				$newpagetemplateTitle = Title::newFromText($new_page_template);
				$templatePage = new WikiPageTE($newpagetemplateTitle);
			
				#gets the content for the template page
				$text1 = $templatePage->getContent();
				$reason = 'Page creation by '.__CLASS__;
			
				#creates a url to link to late
				$t = Title::newFromText($new_page_pageName); 
				$t->getFullURL();
			
				#creates the new page and saves it
				$wikiPage = new WikiPageTE($t);
				$wikiPage->save($text1, $reason);
				
				//$wikiPage->touch(); refreshes the page	
				$wikiPage->touch();	
				
				//Output the created pages to CreatedECKPages.txt. => This file can be used to deleted the crated pages as follows (By using the Mediawiki built-in deleteBatch.php):
				//php /Library/WebServer/Documents/omp/peter/maintenance/deleteBatch.php -u peterwu -r "test" /Users/peterwu/tamu_github/PhpTestCodes/createdECKPages.txt  
				file_put_contents("createdECKPages.txt","OMP_ST:"."$id"."_!_"."$x3\n",FILE_APPEND);
				
				
				
				//new row on strain info table and inserts strain name, parent genotype, and parent name
				$strain_table_template = "Strain_info_table";
				//This pulls out the strain info table.
				$newtable = $wikiPage->getTable($strain_table_template);
				$box = $newtable[0];
				$newrow = $box->insert_row(
					#Strain name
					"Escherichia coli "."$x3"."||".
					#Synonyms
					"||".
					#Taxon Information
					"*Pangenome: Escherichia coli\n*Subspecies and/or strain: ".'K-12'."\n*NCBI Taxonomy ID: ".'[http://www.ncbi.nlm.nih.gov/taxonomy?term=83333 83333]'."||".
					#Genotype
					'?'."||".
					#Strain Reference
					'PMID: 21185072'."||".
					#Strain availability
					"?"."||".
					#Ancestry
					"ancestor:OMP ST:802 ! Escherichia coli K-12 P90C"."||"
					);
				$wikiPage->touch();
				
				//This pulls out OMP_annotation_table
				/*In order for this to work, I might have to enable some other extensions
				$strain_table_template = "OMP_annotation_table";
				$newtable = $new_page->getTable($strain_table_template);
				$box = $newtable[0];
				$newrow = $box->insert_row(
							"|".
							"|".
							"|".
							"|".
							"OMP:(varies)".
							"|".
							"|".
							"abolished viability".
							"|".
							"|".
							"OMP_AN:1028".
							"|".
							"|".
							"|".
							"|".
							"|".
							"|".
							"|".
							"|".
							$reference.
							"|".
							"|"
						);
						$new_page->touch();
				*/
			
			
			}
			else
			{
			echo "$x3 was found!"."\n";
			$j++;
			}
		}
		$k=$j-$i;
		echo "$j strain(s) were found in chemgen under scorelimit = $ScoreLimit, $k out of $j was(were) found to have strain page(s), \n";
		echo "$i page(s) was(were) created\n";
		
		
	
								
			
			
							
		}
				
	//(?)don't understand why this function has to exist
	private function parse_parameters(){
		$this->addOption( "condition", "", $required = false, $withArg = true, $condition = 'c');
		$this->addOption( "strain", "", $required = false, $withArg = true, $strain='s' );
		$this->addOption( "ScoreLimit", "", $required = true, $withArg = true, $ScoreLimit='l' );	//cannot use sc to represent ScoreLimit. Changed $required to true to make sure ScoreLimit has been specified.
		
		/*
		Before I was using this code, but now I undertand that addOption() already does the job:
				
		if($ScoreLimit=="")
		{		
		die ("Must specify Score Limit!\n\n");
		}
		*/
		
		
	}
}

require_once( RUN_MAINTENANCE_IF_MAIN );

/*
Function to set the global variable $IP and include the abstract
class Maintenance
*/
/*
Function to set the global variable $IP and include the abstract
class Maintenance
*/
function set_IP($path){
	global $IP;
	if ( isset($path) && is_file("$path/maintenance/Maintenance.php") ){
		$IP = $path;
		require_once( $IP . "/maintenance/Maintenance.php" );
		return $path; 
	} else {
		//I set my wiki as the default path
		require_once( "/Library/WebServer/Documents/omp/peter/maintenance/Maintenance.php" );
		echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		echo "Default path was used: /Library/WebServer/Documents/omp/peter/maintenance/Maintenance.php\n";
		echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	}
	
}