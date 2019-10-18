<?php

/*
 *  	Sandy Labonte's OMP dev. June 2015
 */

if( defined( 'MW_INSTALL_PATH' ) ) {
	$IP = MW_INSTALL_PATH;
} else {
	$IP = dirname( __FILE__ );
}

$path = array( $IP, "$IP/includes", "$IP/languages" );
set_include_path( implode( PATH_SEPARATOR, $path ) . PATH_SEPARATOR . get_include_path() );

require_once( "$IP/includes/DefaultSettings.php" );

# If PHP's memory limit is very low, some operations may fail.
# ini_set( 'memory_limit', '20M' );

if ( $wgCommandLineMode ) {
	if ( isset( $_SERVER ) && array_key_exists( 'REQUEST_METHOD', $_SERVER ) ) {
		die( "This script must be run from the command line\n" );
	}
}
## Uncomment this to disable output compression
# $wgDisableOutputCompression = true;

$wgSitename         = "OMPwiki";
$wgLogo = '';

## The URL base path to the directory containing the wiki;
## defaults for all runtime URL paths are based off of this.
## For more information on customizing the URLs please see:
## http://www.mediawiki.org/wiki/Manual:Short_URL

# base path from microbialphenotypes.org, net, or com
$wgScriptPath       = "/sandy";
$wgScriptExtension  = ".php";

## UPO means: this is also a user preference option

$wgEnableEmail      = true;
$wgEnableUserEmail  = true; # UPO

$wgEmergencyContact = "jimhu@tamu.edu";
$wgPasswordSender = "jimhu@tamu.edu";
$wgPasswordReminderResendTime = 0; 

$wgEnotifUserTalk = true; # UPO
$wgEnotifWatchlist = true; # UPO
$wgEmailAuthentication = true;

## Database settings
$wgDBtype           = "mysql";
$wgDBserver         = "localhost";
$wgDBname           = "omp_sandy_wikidb";
$wgDBuser           = "wikiuser";
$wgDBpassword       = "dist0phpiati";

# MySQL specific settings
$wgDBprefix         = "";

# MySQL table options to use during installation or update
$wgDBTableOptions   = "ENGINE=InnoDB, DEFAULT CHARSET=binary";

# Experimental charset support for MySQL 4.1/5.0.
$wgDBmysql5 = true;

## Shared memory settings
$wgMainCacheType = CACHE_MEMCACHED;
#$wgMemCachedServers = array( "127.0.0.1:11000" );
$wgMemCachedServers = array( "127.0.0.1:11212" );

## To enable image uploads, make sure the 'images' directory
## is writable, then set this to true:
$wgEnableUploads       = true;
$wgUseImageResize      = true;
# $wgUseImageMagick = true;
# $wgImageMagickConvertCommand = "/usr/bin/convert";

## If you use ImageMagick (or any other shell command) on a
## Linux server, this will need to be set to the name of an
## available UTF-8 locale
# $wgShellLocale = "en_US.UTF-8";

## If you want to use image uploads under safe mode,
## create the directories images/archive, images/thumb and
## images/temp, and make them all writable. Then uncomment
## this, if it's not already uncommented:
# $wgHashedUploadDirectory = false;

## If you have the appropriate support software installed
## you can enable inline LaTeX equations:
$wgUseTeX           = false;

$wgLocalInterwiki   = strtolower( $wgSitename );

$wgLanguageCode = "en";

$wgSecretKey = "4e9549fc80f63e66d8f441ff51804a0f9ed34109f32ae3b53625e4a89b4f3ec6";

## Default skin: you can change the default skin. Use the internal symbolic
## names, ie 'standard', 'nostalgia', 'cologneblue', 'monobook':
$wgDefaultSkin = 'monobook';

## For attaching licensing metadata to pages, and displaying an
## appropriate copyright notice / icon. GNU Free Documentation
## License and Creative Commons licenses are supported so far.
$wgEnableCreativeCommonsRdf = true;
$wgRightsPage = ""; # Set to the title of a wiki page that describes your license/copyright
$wgRightsUrl = "http://www.gnu.org/licenses/old-licenses/fdl-1.2.txt";
$wgRightsText = "GNU Free Documentation License 1.2";
$wgRightsIcon = "${wgScriptPath}/skins/common/images/gnu-fdl.png";
# $wgRightsCode = "gfdl1_2"; # Not yet used

$wgDiff3 = "/usr/bin/diff3";

# When you make changes to this configuration file, this will make
# sure that cached pages are cleared.
$wgCacheEpoch = max( $wgCacheEpoch, gmdate( 'YmdHis', @filemtime( __FILE__ ) ) );

#prevent user account creation

$wgGroupPermissions['*'	   ]['createaccount']   = false;
$wgGroupPermissions['*'    ]['read']            = true;         // 1.5.0
$wgGroupPermissions['*'    ]['edit']            = false;         // 1.5.0
$wgGroupPermissions['*'    ]['createpage']      = false;        // 1.6.0
$wgGroupPermissions['*'    ]['createtalk']      = true;         // 1.6.0
#$wgGroupPermissions['user' ]['createaccount']   = true;
$wgGroupPermissions['user' ]['createpage']   	= true;
$wgGroupPermissions['student']['createaccount']	= false;
$wgGroupPermissions['demo']['edit'] 			= false;
$wgGroupPermissions['demo' ]['createpage']   	= false;
// bot permissions needed for the daily updates
$wgGroupPermissions['bot']['edit']				= true;
$wgGroupPermissions['bot']['createpage']		= true;
$wgGroupPermissions['bot']['createtalk']		= true;
$wgGroupPermissions['bot']['move']				= true;

# Allow linking images
$wgAllowExternalImages = true;

# Allow lowercase page names (need for gene pages)
$wgCapitalLinks =false;

# Add category pages to default search
$wgNamespacesToBeSearchedDefault[NS_CATEGORY] = true;

#allow + in pagename
$wgLegalTitleChars .= '+';

# Allow trackbacks
$wgUseTrackbacks = false;


$wgEmailAuthentication = true;

# https://bugzilla.wikimedia.org/show_bug.cgi?id=11533

# turn off dumb page counters
$wgDisableCounters = true;

# turn off intensive database querying for special pages
$wgMiserMode = false;

#overwrite array of allowed file extensions
$wgFileExtensions = array( 'png', 'tif', 'tiff', 'gif', 'jpg', 'jpeg','pdf', 'txt','pdb', 'doc', 'ppt','pptx', 'svg','pse' );

#$wgReadOnly = 'read only mode for maintenance';
$wgReadOnly = false;

// ======= Code Library ================
require_once( '/usr/local/github/wiki-code/Setup.php' );

// ======= Extensions ==================
/*
	All the extensions for this instance of Mediawiki come from an abstracted place.
	We have all our extension code in a GitLab working directory  in /usr/local/gitlab/wiki-extensions/
	- JH 20140611

	P.S. Everything is listed below in semi-alphabetical order.
*/
$wgExtensionPath 		= 'extensions.local/';
# common external extensions
$wgAbsolutePath2Wiki = dirname(__FILE__);
#require_once("$wgExtensionPath/_extensionsets/commonExternal.php");

$wgAbsolutePath2Wiki = dirname(__FILE__);

$wgImageTmpDir = "/Library/WebServer/Documents/omp/images";
$wgImageTmpURL = "http://microbialphenotypes.org/images";


# Charinsert
    require_once(  $wgExtensionPath . "CharInsert/CharInsert.php" );

# PMID for EUtils
  $extEfetchCache = "/Library/WebServer/tmp/pubmed";
  require_once( $wgExtensionPath . "PMID/PMID.php");

# Cite
	require_once( $wgExtensionPath . "Cite/Cite.php" );
	require_once( $wgExtensionPath . "Cite/ProcessCite.php" );

# Contribution Credits
    if ( $wgCommandLineMode != true ) {
        require_once( $wgExtensionPath . "ContributionCredits_Sidebar.php");
    }
# Datatables
  require_once("extensions/DataTables/DataTables.php");
# ExtendedCategoryPage
    require_once( $wgExtensionPath . 'ExtendedCategoryPage/ExtendedCategoryPage.php' );

# Interwiki
	$wgGroupPermissions['*']['interwiki'] = false;
	$wgGroupPermissions['sysop']['interwiki'] = true;

# OBOfig
    require_once( $wgExtensionPath . "OBOfig.php" );

# OMP
#    require_once( $wgExtensionPath . "/OMP/OMP.php");
#    require_once(  "/usr/local/github/wiki-extensions/wiki-extensions-omp/OMP/OMP.php");

# PagesOnDemand
	require_once( $wgExtensionPath . "PagesOnDemand/PagesOnDemand.php");
	require_once( $wgExtensionPath . "PagesOnDemand/PMID_OnDemand.php");
	require_once( $wgExtensionPath . "PagesOnDemand/GOpageOnDemand.php");
	require_once( $wgExtensionPath . "PagesOnDemand/PreloadOnDemand.php");

# PMID_Summary
    require_once( $wgExtensionPath . "pmid_summary.php"  );


# ProtectSection
	require_once( $wgExtensionPath . "ProtectSection/ProtectSection.php");

# TableEdit
	require_once(  "extensions/TableEdit/TableEdit.php");

# TableEdit modules
	require_once(  "extensions/TableEdit/modules/GONUTS_TableEditLinks.php" );
	require_once( "extensions/TableEdit/modules/TableMarkerUpper.php");
	#require_once(  "extensions/TableEdit/modules/TableEditCreationLinks.php");

# TableEdit column rules
	require_once(  "extensions/TableEdit/column_rules/GO_nr.php");
	require_once(  "extensions/TableEdit/column_rules/strains.php");
	require_once(  "extensions/TableEdit/column_rules/column_rules_settings.php");
	#require_once( "/Volumes/home/shabnam/wiki-extensions_trunk/TableEdit/column_rules/column_rules_settings.php");

	require_once( "extensions/TableEdit/column_rules/dbxref.php" );
# UserRightsList
	require_once( $wgExtensionPath . "UserRightsList/UserRightsList.php");

# User Recent
	require_once( $wgExtensionPath . "UserRecent.php");

# Restrictive Rights
	require_once( $wgExtensionPath . "RestrictiveRights.php");

# RSS News
	require_once( $wgExtensionPath . "RSSNews.php");

# OMP summary
	require_once( $wgExtensionPath . "OMP/OMPsummary.php");

#include("extension/GraphViz/GraphViz.php"); //includes the extension
$wgGraphVizSettings->execPath = "/path/to/graphviz/"; //Path to the Graphviz Installation

/**
 *  Daniel wrote this next function to format the categories at the
 *  bottom of the page in a more readable way.

if ($wgVersion < "1.19"){
   $wgHooks['OutputPageMakeCategoryLinks'][] = 'efBetterCategoriesSectionInArticles';
}
function efBetterCategoriesSectionInArticles( &$out, $categories, &$mCategoryLinks ) {
	global $wgContLang;


	$categories = array_unique( array_keys($categories) );

	if ( $categories ) {
	
		function compareGoTerms( $a, $b ) {
			// GO:0004713_!_
			// 1234567890123
			if ( preg_match( '/^GO:\d+_!_(.*)/', $a, $m ) ) {
				$a = substr($m[1], 0, 12);
			}
			if ( preg_match( '/^GO:\d+_!_(.*)/', $b, $m ) ) {
				$b = substr($m[1], 0, 12);
			}
			return strnatcasecmp( $a, $b );
		}
	
		$c = new CategoryViewer( Title::newFromText( $categories[0], NS_CATEGORY ) );
		$l = new Linker;
		usort($categories, 'compareGoTerms');
		$page = array();
		$page_start_char = array();
		foreach ( $categories as $category )  {
			$sortkey = strtolower( $category );	
			// ..for GO terms
			if ( preg_match( '/^GO:\d+_!_(.*)/', $category, $m ) ) {
				$sortkey = strtolower( $m[1] );
			}
			$t = Title::newFromText( $category, NS_CATEGORY );
			$page[] = '<span class="redirect-in-category">' . $l->link($t, $t->getText()) . '</span>';
			$page_start_char[] = $wgContLang->convert( $wgContLang->firstChar( $sortkey ) );
		}		
		usort($page_start_char,  'strnatcasecmp');
		$html = $c->formatList( $page, $page_start_char );			
		
		 // 	
		 //	Here you could really format it however you wanted...with
		 //	  some nice JavaScript or a jQuery library.
		 //
		 //	$html = '<div class="betterCategories">' . . '</div>';
 		 //   ...etc...
		
	}

	$html = '<div style="font-size:smaller;">' .  $html . '</div>';

	$mCategoryLinks['normal'] = array( $html );

	return false;
}
*/





// ======== rotate $wgLogo =========================================




#Block msnbot
$agent= " " . @$_SERVER['HTTP_USER_AGENT'];
if ( strpos($agent,"msnbot") ) {
        exit;
}
$wgShowExceptionDetails = true;
$wgShowSQLErrors = true;

# requires for Sandy's customization
#require_once('/Volumes/home/sandylabonte');
#require_once("/Volumes/home/sandylabonte/phptutorial/Unit_3/SampleExtension/SampleExtension.php");
#calling samplespecialpage
#require_once("/Volumes/home/sandylabonte/phptutorial/Unit_3/SampleSpecialPage/SampleSpecialPage.php" );
#require_once("/Users/sandylabonte/github/wiki-extensions-omp/OMP/OMP.php" );
require_once("/Users/beckyberg/github/wiki-extensions-omp/OMP/BeckySpecialpage/beckysamplespecialpage.php");
require_once("/Users/beckyberg/github/wiki-extensions-omp/OMP/EcoliSpecialpage/Ecolisamplespecialpage.php");
require_once("/Users/sandylabonte/github/wiki-extensions-omp/OMP/ECGASpecialpage/ECGAsamplespecialpage.php");
require_once("/usr/local/github/wiki-extensions/wiki-extensions-mw/SimpleTooltip/SimpleTooltip.php");