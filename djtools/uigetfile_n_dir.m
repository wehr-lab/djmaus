<<<<<<< HEAD
function [pathname] = uigetdir2(start_path, dialog_title)
=======
function [pathname] = uigetfile_n_dir(start_path, dialog_title)
>>>>>>> f2143cfd74a66f9f3877bb57c8a2c98ede027402
% Pick multiple directories and/or files

import javax.swing.JFileChooser;

<<<<<<< HEAD
=======
% if nargin == 0 || start_path == '' || start_path == 0 % Allow a null argument.
%     start_path = pwd;
% end
>>>>>>> f2143cfd74a66f9f3877bb57c8a2c98ede027402

jchooser = javaObjectEDT('javax.swing.JFileChooser', start_path);

jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
if nargin > 1
    jchooser.setDialogTitle(dialog_title);
end

jchooser.setMultiSelectionEnabled(true);

status = jchooser.showOpenDialog([]);

if status == JFileChooser.APPROVE_OPTION
    jFile = jchooser.getSelectedFiles();
	pathname{size(jFile, 1)}=[];
    for i=1:size(jFile, 1)
		pathname{i} = char(jFile(i).getAbsolutePath);
	end
	
elseif status == JFileChooser.CANCEL_OPTION
    pathname = [];
else
    error('Error occured while picking file.');
end
