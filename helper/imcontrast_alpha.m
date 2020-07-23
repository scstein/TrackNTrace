function him = imcontrast_alpha(ax)
% him = IMCONTRAST_ALPHA(ax) wraper to use imcontrast(ax) for ALim instead of CLim
% Creates a invisible axes with the alpha data on top of the given axes and synchronises the CLim with ALim.
%
% Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de, 2020

	if nargin==0 || isempty(ax)
		ax = gca;
	end
	
	[imageHandle, axHandle, figHandle] = imhandles(ax);
	imageHandle = imageHandle(1); %Take the first image in case of multiple images
	% make an invisble copy of the axes and plot the alpha channel
	ax_helper = axes(figHandle,'Position',get(axHandle,'Position'),'Visible','off');
%         uistack(ax_helper,'bottom') 
	
	imagehandle_helper = imagesc(ax_helper,get(imageHandle,'AlphaData'),'AlphaData',0);
	set(ax_helper ,'Visible','off')
	set(ax_helper,'CLim',get(axHandle,'ALim'));
	set(imagehandle_helper,'CDataMapping',get(imageHandle,'AlphaDataMapping'));
	
	limlistener = addlistener(ax_helper,'CLim','PostSet',@(src,evnt)set(axHandle,'ALim',get(ax_helper,'CLim')));
	datalistener = addlistener(imagehandle_helper,'CData','PostSet',@(src,evnt)set(imageHandle,'AlphaData',get(imagehandle_helper,'CData')));
	
	him = imcontrast(ax_helper);
	set(him,'CloseRequestFcn',@cleanUp)
	
	function cleanUp(~,~)
		try
			delete(limlistener);
			delete(datalistener);
			delete(ax_helper);
		catch err
			warning(err.message);
		end
	end
end