function fullPath = getFullPath( relPath )
curPath = pwd() ;
cd(relPath) ;
fullPath = pwd() ;
cd( curPath ) ;
