classdef app_conservation_laws_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        GridLayout           matlab.ui.container.GridLayout
        LeftPanel            matlab.ui.container.Panel
        GridLayout3          matlab.ui.container.GridLayout
        FluxListBox          matlab.ui.control.ListBox
        FluxListBoxLabel     matlab.ui.control.Label
        uREditField          matlab.ui.control.NumericEditField
        uREditFieldLabel     matlab.ui.control.Label
        uLEditField          matlab.ui.control.NumericEditField
        uLEditFieldLabel     matlab.ui.control.Label
        dxEditField          matlab.ui.control.NumericEditField
        dxEditFieldLabel     matlab.ui.control.Label
        MetodoListBox        matlab.ui.control.ListBox
        MetodoListBoxLabel   matlab.ui.control.Label
        TimeEditField        matlab.ui.control.NumericEditField
        TimeEditFieldLabel   matlab.ui.control.Label
        StartButton          matlab.ui.control.Button
        RightPanel           matlab.ui.container.Panel
        GridLayout2          matlab.ui.container.GridLayout
        dtEditField          matlab.ui.control.EditField
        dtEditFieldLabel     matlab.ui.control.Label
        RiemannProblemLabel  matlab.ui.control.Label
        UIAxes               matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = private)
        T = 1;
        dx = 0.01
        uL = 0.2;
        uR = 0.9;
        id = 1;
        fb = 'Burgers';
        Us;
        x;
        dt;
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)

            app.UIFigure.Position = [200, 200, 768, 576];
            movegui(app.UIFigure, 'center');
            %axis(app.UIAxes, 'square');
            app.TimeEditField.Value = app.T;
            app.dxEditField.Value = app.dx;
            app.uLEditField.Value = app.uL;
            app.uREditField.Value = app.uR;
        end

        % Value changed function: FluxListBox
        function FluxListBoxValueChanged(app, event)
            app.fb = app.FluxListBox.Value;
        end

        % Callback function
        function StartButtonValueChanged(app, event)
            value = app.StartButton.Value;
            
        end

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
            [app.Us,app.x,app.dt] = riemann_app_script(app.T,app.dx,app.uL,app.uR,app.fb,app.id);
            app.dtEditField.Value = string(app.dt);
            p = plot(app.UIAxes,app.x,app.Us(:,1));
            legend(app.UIAxes,"u solution");
            title(app.UIAxes,"T = " + 0)
            ylim(app.UIAxes,[min(app.uL,app.uR) - 0.2,max(app.uL,app.uR) + 0.2])
            for n = 2:size(app.Us,2)
                p.YData = app.Us(:,n);
                title(app.UIAxes,"T = " + (n-1)*app.dt)
                drawnow
                pause(0.03)
            end
        end

        % Value changed function: TimeEditField
        function TimeEditFieldValueChanged(app, event)

            app.T = app.TimeEditField.Value;
            
        end

        % Value changed function: dxEditField
        function dxEditFieldValueChanged(app, event)
            app.dx = app.dxEditField.Value;
            
        end

        % Value changed function: uREditField
        function uREditFieldValueChanged(app, event)
            app.uR = app.uREditField.Value;
            
        end

        % Value changed function: uLEditField
        function uLEditFieldValueChanged(app, event)
            app.uL = app.uLEditField.Value;
            
        end

        % Clicked callback: MetodoListBox
        function MetodoListBoxClicked(app, event)
            app.id = event.InteractionInformation.Item;
            
        end

        % Callback function
        function TempofinaleEditFieldValueChanged2(app, event)
            
           
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {480, 480};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {224, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Theme = 'light';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);
            app.UIFigure.WindowStyle = 'modal';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {224, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.LeftPanel);
            app.GridLayout3.ColumnWidth = {37, 26, 88};
            app.GridLayout3.RowHeight = {22, '1.37x', 22, 22, 22, 22, 22, '3.24x', '1x', 23};

            % Create StartButton
            app.StartButton = uibutton(app.GridLayout3, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.FontSize = 24;
            app.StartButton.FontWeight = 'bold';
            app.StartButton.Layout.Row = 9;
            app.StartButton.Layout.Column = [2 3];
            app.StartButton.Text = 'Start';

            % Create TimeEditFieldLabel
            app.TimeEditFieldLabel = uilabel(app.GridLayout3);
            app.TimeEditFieldLabel.HorizontalAlignment = 'center';
            app.TimeEditFieldLabel.FontSize = 14;
            app.TimeEditFieldLabel.Layout.Row = 3;
            app.TimeEditFieldLabel.Layout.Column = 1;
            app.TimeEditFieldLabel.Text = 'Time';

            % Create TimeEditField
            app.TimeEditField = uieditfield(app.GridLayout3, 'numeric');
            app.TimeEditField.ValueChangedFcn = createCallbackFcn(app, @TimeEditFieldValueChanged, true);
            app.TimeEditField.HorizontalAlignment = 'center';
            app.TimeEditField.FontSize = 14;
            app.TimeEditField.Layout.Row = 3;
            app.TimeEditField.Layout.Column = [2 3];

            % Create MetodoListBoxLabel
            app.MetodoListBoxLabel = uilabel(app.GridLayout3);
            app.MetodoListBoxLabel.HorizontalAlignment = 'center';
            app.MetodoListBoxLabel.Layout.Row = 7;
            app.MetodoListBoxLabel.Layout.Column = 1;
            app.MetodoListBoxLabel.Text = 'Metodo';

            % Create MetodoListBox
            app.MetodoListBox = uilistbox(app.GridLayout3);
            app.MetodoListBox.Items = {'Upwind', 'Lax-F', 'MacCormack', 'Lax-W', 'Upwind correction convex case', 'Upwind correction concave case', 'Godunov convex', 'Godunov concave'};
            app.MetodoListBox.FontSize = 14;
            app.MetodoListBox.Layout.Row = [7 8];
            app.MetodoListBox.Layout.Column = [2 3];
            app.MetodoListBox.ClickedFcn = createCallbackFcn(app, @MetodoListBoxClicked, true);
            app.MetodoListBox.Value = 'Upwind';

            % Create dxEditFieldLabel
            app.dxEditFieldLabel = uilabel(app.GridLayout3);
            app.dxEditFieldLabel.HorizontalAlignment = 'center';
            app.dxEditFieldLabel.FontSize = 14;
            app.dxEditFieldLabel.Layout.Row = 4;
            app.dxEditFieldLabel.Layout.Column = 1;
            app.dxEditFieldLabel.Text = 'dx';

            % Create dxEditField
            app.dxEditField = uieditfield(app.GridLayout3, 'numeric');
            app.dxEditField.ValueChangedFcn = createCallbackFcn(app, @dxEditFieldValueChanged, true);
            app.dxEditField.HorizontalAlignment = 'center';
            app.dxEditField.FontSize = 14;
            app.dxEditField.Layout.Row = 4;
            app.dxEditField.Layout.Column = [2 3];

            % Create uLEditFieldLabel
            app.uLEditFieldLabel = uilabel(app.GridLayout3);
            app.uLEditFieldLabel.HorizontalAlignment = 'center';
            app.uLEditFieldLabel.FontSize = 14;
            app.uLEditFieldLabel.Layout.Row = 6;
            app.uLEditFieldLabel.Layout.Column = 1;
            app.uLEditFieldLabel.Text = 'uL';

            % Create uLEditField
            app.uLEditField = uieditfield(app.GridLayout3, 'numeric');
            app.uLEditField.ValueChangedFcn = createCallbackFcn(app, @uLEditFieldValueChanged, true);
            app.uLEditField.HorizontalAlignment = 'center';
            app.uLEditField.FontSize = 14;
            app.uLEditField.Layout.Row = 6;
            app.uLEditField.Layout.Column = [2 3];

            % Create uREditFieldLabel
            app.uREditFieldLabel = uilabel(app.GridLayout3);
            app.uREditFieldLabel.HorizontalAlignment = 'center';
            app.uREditFieldLabel.FontSize = 14;
            app.uREditFieldLabel.Layout.Row = 5;
            app.uREditFieldLabel.Layout.Column = 1;
            app.uREditFieldLabel.Text = 'uR';

            % Create uREditField
            app.uREditField = uieditfield(app.GridLayout3, 'numeric');
            app.uREditField.ValueChangedFcn = createCallbackFcn(app, @uREditFieldValueChanged, true);
            app.uREditField.HorizontalAlignment = 'center';
            app.uREditField.FontSize = 14;
            app.uREditField.Layout.Row = 5;
            app.uREditField.Layout.Column = [2 3];

            % Create FluxListBoxLabel
            app.FluxListBoxLabel = uilabel(app.GridLayout3);
            app.FluxListBoxLabel.HorizontalAlignment = 'center';
            app.FluxListBoxLabel.FontSize = 14;
            app.FluxListBoxLabel.Layout.Row = 1;
            app.FluxListBoxLabel.Layout.Column = 1;
            app.FluxListBoxLabel.Text = 'Flux';

            % Create FluxListBox
            app.FluxListBox = uilistbox(app.GridLayout3);
            app.FluxListBox.Items = {'Burgers', 'Traffic'};
            app.FluxListBox.ValueChangedFcn = createCallbackFcn(app, @FluxListBoxValueChanged, true);
            app.FluxListBox.FontSize = 14;
            app.FluxListBox.Layout.Row = [1 2];
            app.FluxListBox.Layout.Column = [2 3];
            app.FluxListBox.Value = 'Burgers';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.RightPanel);
            app.GridLayout2.ColumnWidth = {'1x', 210, '1x'};
            app.GridLayout2.RowHeight = {31, '1.03x', '3.49x', '1x'};
            app.GridLayout2.ColumnSpacing = 2.5;
            app.GridLayout2.Padding = [2.5 10 2.5 10];

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout2);
            title(app.UIAxes, 'Evolution')
            xlabel(app.UIAxes, 'x')
            ylabel(app.UIAxes, 'u')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontName = 'Hack';
            app.UIAxes.Layout.Row = 3;
            app.UIAxes.Layout.Column = [1 3];

            % Create RiemannProblemLabel
            app.RiemannProblemLabel = uilabel(app.GridLayout2);
            app.RiemannProblemLabel.HorizontalAlignment = 'center';
            app.RiemannProblemLabel.FontName = 'IBM Plex Mono';
            app.RiemannProblemLabel.FontSize = 24;
            app.RiemannProblemLabel.FontWeight = 'bold';
            app.RiemannProblemLabel.FontColor = [0.0667 0.4431 0.7451];
            app.RiemannProblemLabel.Layout.Row = 1;
            app.RiemannProblemLabel.Layout.Column = [1 3];
            app.RiemannProblemLabel.Text = 'Riemann Problem';

            % Create dtEditFieldLabel
            app.dtEditFieldLabel = uilabel(app.GridLayout2);
            app.dtEditFieldLabel.HorizontalAlignment = 'center';
            app.dtEditFieldLabel.FontSize = 18;
            app.dtEditFieldLabel.Layout.Row = 4;
            app.dtEditFieldLabel.Layout.Column = 1;
            app.dtEditFieldLabel.Text = 'dt';

            % Create dtEditField
            app.dtEditField = uieditfield(app.GridLayout2, 'text');
            app.dtEditField.HorizontalAlignment = 'center';
            app.dtEditField.FontSize = 18;
            app.dtEditField.Layout.Row = 4;
            app.dtEditField.Layout.Column = 2;
            app.dtEditField.Value = '0';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app_conservation_laws_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end