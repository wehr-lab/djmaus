function age_out=age(dob, day)
%returns age of animal in days
%usage: 
%     age(dob) to get animal's age today
%     age(dob, date) to get animal's age on a particular date

if nargin==1
    age_out=datenum(date)-datenum(dob);
elseif nargin==2
    age_out=round(datenum(day)-datenum(dob));
end