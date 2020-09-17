function DateNow = bgc1d_getDate()
 DateNow = datestr(now);
 idxreplace=[strfind(DateNow,':'),strfind(DateNow,' '),strfind(DateNow,'-')];
 DateNow(idxreplace)='_';
